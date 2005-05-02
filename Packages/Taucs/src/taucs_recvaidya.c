/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*                                                       */
/* Recursive Vaidya preconditioners.                     */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

/*#include <unistd.h>*/

/*
long int random();
void srandom(unsigned int seed);
*/

#ifdef TAUCS_CORE

typedef struct {
  taucs_ccs_matrix** B;
  taucs_ccs_matrix** S;
  taucs_ccs_matrix** L;
  int             levels;
  int             level;
  double          convratio;
  double          maxits;
} recvaidya_args;

static int
recvaidya_order(taucs_ccs_matrix* m,
		int** perm,
		int** invperm,
		int*  P)
{
  int  n,nnz,i,j,ip,k,p,nleaves;
  int* adjptr;
  int* adj;
  int* len;
  int* ptr;
  int* degree;
  int* leaves;

  n   = m->n;
  nnz = (m->colptr)[n];
  
  taucs_printf("recvaidya_order: starting, matrix is %dx%d, # edges=%d\n",
	     n,n,nnz-n);

  *perm    = (int*) taucs_malloc(n * sizeof(int));
  *invperm = (int*) taucs_malloc(n * sizeof(int));

  /* we can reuse buffers: don't need invperm until the end */
  /* also, we can reuse perm for leaves but it's messy.     */
  len    = (int*) taucs_malloc(n * sizeof(int));
  degree = (int*) taucs_malloc(n * sizeof(int));
  leaves = (int*) taucs_malloc(n * sizeof(int));

  adjptr = (int*) taucs_malloc(n * sizeof(int));
  adj    = (int*) taucs_malloc(2*(nnz-n) * sizeof(int));

  if (!(*perm) || !(*invperm) || !adjptr || !adj || !len || !degree || ! leaves) {
    taucs_free(adj);
    taucs_free(adjptr);
    taucs_free(len);
    taucs_free(leaves);
    taucs_free(degree);
    taucs_free(*perm);
    taucs_free(*invperm);
    return -1;
  }

  for (i=0; i<n; i++) len[i] = 0;

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      /*i = (m->rowind)[ip] - (m->indshift);*/
      i = (m->rowind)[ip];
      if (i != j) {
	len[i] ++;
	len[j] ++;
      }
    }
  }

  nleaves = 0;
  for (i=0; i<n; i++) {
    degree[i] = len[i]; 
    if (degree[i] <= 1) {
      leaves[nleaves] = i;
      nleaves++;
    }
  }

  adjptr[0] = 0;
  for (i=1; i<n; i++) adjptr[i] = adjptr[i-1] + len[i-1];

  ptr =  *perm;
  for (i=0; i<n; i++) ptr[i] = adjptr[i];

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      /*i = (m->rowind)[ip] - (m->indshift);*/
      i = (m->rowind)[ip];
      if (i != j) {
	adj[ ptr[i] ] = j;
	adj[ ptr[j] ] = i;
	ptr[i] ++;
	ptr[j] ++;
      }
    }
  }

  /* now the graph data structure is ready */

  /* we first eliminate leaves until all the degrees >= 2 */

  i = 0;
  while (nleaves > 0) {
    /*taucs_printf("recvaidya_order: nleaves=%d\n",nleaves);*/
    nleaves--;
    j = leaves[nleaves];

    /*taucs_printf("recvaidya_order: next leaf is %d, degree=%d\n",j,len[j]);*/
    
    (*perm)   [ i ] = j;
    (*invperm)[ j ] = i;
    i++;
    
    if (len[j] > 0) {
      assert(len[j] == 1); /* j must be a degree-1 vertex */
      len[j]--;            /* eliminate j */
      p = adj[ adjptr[j] ]; 
      
      /*taucs_printf("symccs_treeorder: parent of %d is %d\n",j,p);*/

      for (k = 0; k < len[p]; k++)
	if (adj[ adjptr[p] + k ] == j) break;

      assert( k < len[p] ); /* j must be a neighbor of p */
	
      /* now delete j from p's adjacency list and compress */
      len[p] --;
      for (; k < len[p]; k++)
	adj[ adjptr[p] + k ] = adj[ adjptr[p] + k+1 ];

      if (len[p] == 1) { /* degree was higher and now is 1 */
	leaves[ nleaves ] = p;
	nleaves++;
      }
    }
  }

  /* an eliminated vertix j must have len[j]==0        */
  /* we can now eliminate all the degree-2 vertices    */
  /* elimination of degree-2 vertices does not change  */
  /* degrees, so we first find them and then eliminate */

  for (j=0; j<n; j++) {
    if (len[j] == 2) {
      (*perm)[i]    = j;
      (*invperm)[j] = i;
      i++;
      len[j] = 0; /* eliminate from the graph */
    }
  }

  *P = i;
  taucs_printf("recvaidya_order: eliminating %d vertices (remaining have deg>2)\n",
	     *P);
  
  for (j=0; j<n; j++) {
    if (len[j] > 0) {
      (*perm)[i]    = j;
      (*invperm)[j] = i;
      i++;
    }
  }

  assert( i == n );

  taucs_free(adj);
  taucs_free(adjptr);
  taucs_free(len);
  taucs_free(leaves);
  taucs_free(degree);

  taucs_printf("recvaidya_order: done\n");

  return 0;
}

void* 
taucs_recursive_amwb_preconditioner_create(taucs_ccs_matrix* A, 
					   double c, 
					   double epsilon, 
					   int nsmall,
					   int maxlevels,
					   int innerits,
					   double convratio,
					   int** perm, 
					   int** invperm)
{
  int l,i,k;
  int levels;
  int P[32];

#if 0
  taucs_ccs_matrix* Sx[32]; /* Schur complements                        */
  taucs_ccs_matrix* Lx[32]; /* Partial LL^T factors                     */
#endif

  taucs_ccs_matrix** S; /* Schur complements                        */
  taucs_ccs_matrix** L; /* Partial LL^T factors                     */

  double exponent = 1.0/(1.0+epsilon);

  recvaidya_args* args;

  int* perml;  /* local permutation for level l */
  int* iperml; /* local permutation for level l */
  int* tmpperm;

  int next = 0;

  if (maxlevels > 32) {
    taucs_printf("taucs_recursive_amwb_preconditioner_create: maxlevel must be 32 or less\n");
    return NULL;
  }

  args = (recvaidya_args*) taucs_malloc(sizeof(recvaidya_args));
  S    = (taucs_ccs_matrix**) taucs_malloc(32 * sizeof(taucs_ccs_matrix*));
  L    = (taucs_ccs_matrix**) taucs_malloc(32 * sizeof(taucs_ccs_matrix*));

  *perm    = (int*) taucs_malloc(A->n * sizeof(int));
  *invperm = (int*) taucs_malloc(A->n * sizeof(int));
  tmpperm = *invperm;
  assert(args && *perm && *invperm);

  for (i=0; i<A->n; i++) (*perm)[i] = (*invperm)[i] = i;

  for (l=0; l<32; l++)
    S[l] = L[l] = NULL;
  
  for (l=0; l<maxlevels; l++) {

    taucs_ccs_matrix* Al;
    taucs_ccs_matrix* V;    /* a Vaidya preconditioner */
    taucs_ccs_matrix* PVPT;
    taucs_ccs_matrix* Ll;
    int     p;
    int     n;
    int rnd;
    int seed = 123;
    double t;

    if (l==0) Al = A;
    else      Al = S[l];

    n = Al->n;
    if (n==0) {l--; break;}

    t = c * pow( (double)n, exponent );
    taucs_printf("recvaidya_create: n=%d c=%.2le eps=%.2le ==> t=%.0lf\n",
	       n,c,epsilon,t);
    srand(seed);
    rnd = rand();
    V = taucs_amwb_preconditioner_create(Al,rnd,t,0 /* stretch flag */);

   if (n <= nsmall || l==maxlevels-1) {
      taucs_printf("recvaidya_create: n=%d <= nsmall=%d (or max level)\n",
		 n,nsmall);
      taucs_ccs_order(V,&perml,&iperml,"metis");
      PVPT = taucs_ccs_permute_symmetrically(V,perml,iperml);
      taucs_ccs_free(V);
      /*
      taucs_ccs_order(Al,&perml,&iperml,"metis");
      PVPT = taucs_ccs_permute_symmetrically(Al,perml,iperml);
      */
      p=n;
    } else {
      recvaidya_order(V,&perml,&iperml,&p);
      PVPT = taucs_ccs_permute_symmetrically(V,perml,iperml);
      taucs_ccs_free(V);
      /*
      taucs_ccs_order(V,&perml,&iperml,"md");
      */
      if (p>n) p=n;
    }

    P[l] = p;

    /* now compose the permutations */

    for (i=0; i<next; i++)
      tmpperm[i] = (*perm)[i];
    for (i=next; i<next+n; i++)
      tmpperm[i] = (*perm)[ next + perml[ i - next ] ];

    for (i=0; i<next+n; i++)
      (*perm)[i] = tmpperm[i];
    for (i=0; i<next+n; i++)
      (*invperm)[ (*perm)[i] ] = i;

    for (k=1; k<=l; k++) {
      int* backperm;
      int* ibackperm;

      taucs_ccs_matrix* PLPT;
      taucs_ccs_matrix* PSPT;

      int N;

      N = L[k-1]->n; 

      backperm  = (int*) taucs_malloc(N * sizeof(int));
      ibackperm = (int*) taucs_malloc(N * sizeof(int));

      for (i=0   ; i<N-n ; i++) backperm[i] = i;
      for (i=N-n ; i<N   ; i++) backperm[i] = (N-n) + perml[i-(N-n)];
      for (i=0   ; i<N   ; i++) ibackperm[backperm[i]] = i;
      
      PLPT = taucs_ccs_permute_symmetrically(L[k-1],backperm,ibackperm);
      taucs_ccs_free(L[k-1]);
      L[k-1] = PLPT;

      taucs_free(backperm);
      taucs_free(ibackperm);

      N = S[k]->n; 

      backperm  = (int*) taucs_malloc(N * sizeof(int));
      ibackperm = (int*) taucs_malloc(N * sizeof(int));

      for (i=0   ; i<N-n ; i++) backperm[i] = i;
      for (i=N-n ; i<N   ; i++) backperm[i] = (N-n) + perml[i-(N-n)];
      for (i=0   ; i<N   ; i++) ibackperm[backperm[i]] = i;
      
      PSPT = taucs_ccs_permute_symmetrically(S[k],backperm,ibackperm);
      taucs_ccs_free(S[k]);
      S[k] = PSPT;

      taucs_free(backperm);
      taucs_free(ibackperm);
    }

    if (p<n) {
      Ll = taucs_ccs_factor_llt_partial(PVPT,p);
      taucs_ccs_free(PVPT);
      taucs_ccs_split(Ll,&(L[l]),&(S[l+1]),p);
      (L[l])   -> flags = TAUCS_TRIANGULAR | TAUCS_LOWER;
      (S[l+1]) -> flags = TAUCS_SYMMETRIC  | TAUCS_LOWER;
      taucs_ccs_free(Ll);
    } else {
      L[l] = taucs_ccs_factor_llt(PVPT,0.0,0);
      taucs_ccs_free(PVPT);
      break;
    }

    next += p;
  }

  levels = l+1;

  taucs_printf("recvaidya-create: %d levels [ ",levels);
  for (l=0; l<levels; l++) taucs_printf("%d ",P[l]);
  taucs_printf("]\n");

  args->levels = levels;
  args->convratio = convratio;
  args->maxits = innerits;
  args->level  = 0;
  args->S = S;
  args->L = L;

  return args;
  
}

int
taucs_recursive_amwb_preconditioner_solve(void* vP,
					  void* vZ, 
					  void* vR)
{
  recvaidya_args* P = (recvaidya_args*) vP;
  double* Z = (double*) vZ;
  double* R = (double*) vR;
  recvaidya_args args;

  if ( P->level == (P->levels)-1 ) {
    /* this is the last level, L is a complete factor */

    /*
    taucs_printf("recvaidya_solve: level=%d/%d, direct solve\n",
	       P->level,P->levels);
    */
    
    taucs_ccs_solve_llt((P->L)[P->level],
			   Z, R);
  } else {
    /*
    taucs_printf("recvaidya_solve: level=%d/%d, Schur complement solve\n",
	       P->level,P->levels);
    */

    args       = *P; /* copy the data but modify next level! */
    args.level = (P->level) + 1;
    
    taucs_ccs_solve_schur((P->L)[P->level],
			  (P->S)[(P->level) + 1],
			  taucs_recursive_amwb_preconditioner_solve,
			  &args,
			  (int)(P->maxits),
			  (P->convratio),
			  Z, R);
  }
  return 0;
}

#endif /* TAUCS_CORE */  
		      




