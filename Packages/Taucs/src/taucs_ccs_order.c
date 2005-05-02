/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

/*********************************************************/
/* Interface to AMD                                      */
/*********************************************************/

#include "../external/src/colamd.h"

static void
taucs_ccs_colamd(taucs_ccs_matrix* m, 
		 int** perm, int** invperm,
		 char* which)
{
#ifndef TAUCS_CONFIG_COLAMD
  taucs_printf("taucs_ccs_colamd: COLAMD routines not linked.\n");
  *perm    = NULL;
  *invperm = NULL;
  return;
#else
  double knobs[COLAMD_KNOBS];
  int    Alen;
  int*   A;
  int*   p;
  int*   ip;
  int    k,nnz;
  int i;
  
  if (m->flags & TAUCS_SYMMETRIC || m->flags & TAUCS_HERMITIAN) {
    taucs_printf("taucs_ccs_colamd: not applicable for symmetric or hermitian matrices\n");
    return;
  }

  taucs_printf("taucs_ccs_colamd: starting\n");
  
  *perm    = NULL;
  *invperm = NULL;

  nnz = (m->colptr)[m->n];
  
  p  = (int*) taucs_malloc((m->n + 1) * sizeof(int));
  ip = (int*) taucs_malloc((m->n + 1) * sizeof(int));
  assert(p && ip);
  
  Alen = colamd_recommended(nnz, m->m, m->n);
  A = taucs_malloc(Alen * sizeof(int)); assert(A);
  assert(A);
  colamd_set_defaults (knobs) ;
  
  for (i=0; i<=m->n; i++)  p[i] = (m->colptr)[i];
  for (k=0; k<nnz; k++)     A[k] = (m->rowind)[k];
  
  taucs_printf("oocsp_ccs_colamd: calling colamd matrix is %dx%d, nnz=%d\n",
	     m->m,m->n,nnz);
  if (!colamd (m->m, m->n, Alen, A, p, knobs)) {
    taucs_printf("oocsp_ccs_colamd: colamd failed\n");
    taucs_free(A);
    taucs_free(p);
    return;
  }
  taucs_printf("oocsp_ccs_colamd: colamd returned\n");
  
  taucs_free(A);

  *perm    = p;
  *invperm = ip;

  for (i=0; i<m->n; i++) (*invperm)[(*perm)[i]] = i;
#endif
}

/*********************************************************/
/* Interface to AMD                                      */
/*********************************************************/

static void
taucs_ccs_amd(taucs_ccs_matrix* m, 
	      int** perm, int** invperm,
	      char* which)
{
#ifndef TAUCS_CONFIG_AMD
  taucs_printf("taucs_ccs_amd: AMD routines not linked.\n");
  *perm    = NULL;
  *invperm = NULL;
  return;
#else
  int  n, iwlen, pfree, ncmpa, iovflo;
  int* iw;
  int* pe;
  int* degree;
  int* nv;
  int* next;
  int* last;
  int* head;
  int* elen;
  int* w;
  int* len;

  int  nnz,i,j,ip;
  
  taucs_printf("taucs_ccs_amd: starting (%s)\n",which);

  if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
    taucs_printf("taucs_ccs_amd: AMD ordering only works on symmetric matrices.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }
  /* this routine may actually work on UPPER as well */
  if (!(m->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_amd: the lower part of the matrix must be represented.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }
    
  *perm    = NULL;
  *invperm = NULL;

  n   = m->n;
  nnz = (m->colptr)[n];
  
  pe     = (int*) taucs_malloc(n * sizeof(int));
  degree = (int*) taucs_malloc(n * sizeof(int));
  nv     = (int*) taucs_malloc(n * sizeof(int));
  next   = (int*) taucs_malloc(n * sizeof(int));
  last   = (int*) taucs_malloc(n * sizeof(int));
  head   = (int*) taucs_malloc(n * sizeof(int));
  elen   = (int*) taucs_malloc(n * sizeof(int));
  w      = (int*) taucs_malloc(n * sizeof(int));
  len    = (int*) taucs_malloc(n * sizeof(int));

  /* AMD docs recommend iwlen >= 1.2 nnz, but this leads to compressions */
  iwlen = n + (int) (2.0 * 2.0*(nnz - n));

  taucs_printf("taucs_ccs_amd: allocating %d ints for iw\n",iwlen);

  iw = (int*) taucs_malloc(iwlen * sizeof(int));

  if (!pe || !degree || !nv || !next || !last || !head 
      || !elen || !w || !len || !iw) {
    taucs_printf("taucs_ccs_amd: out of memory\n");
    taucs_free(pe    );
    taucs_free(degree);
    taucs_free(nv    );
    taucs_free(next  );
    taucs_free(last  );
    taucs_free(head  );
    taucs_free(elen  );
    taucs_free(w     );
    taucs_free(len   );
    taucs_free(iw    );
    return;
  }

  /*
  assert(iw && pe && degree && nv && next && last && head &&
	 elen && w && len); 
  */

  assert(sizeof(int) == 4);
  /*iovflo = 2147483648; */ /* for 32-bit only! */
  iovflo = 2147483647; /* for 32-bit only! */

  for (i=0; i<n; i++) len[i] = 0;

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      i = (m->rowind)[ip];
      /*i = (m->rowind)[ip] - (m->indshift);*/
      if (i != j) {
	len[i] ++;
	len[j] ++;
      }
    }
  }

  pe[0] = 1;
  for (i=1; i<n; i++) pe[i] = pe[i-1] + len[i-1];
  
  pfree = pe[n-1] + len[n-1];

  /* use degree as a temporary */

  for (i=0; i<n; i++) degree[i] = pe[i] - 1;

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      /*i = (m->rowind)[ip] - (m->indshift);*/
      i = (m->rowind)[ip];
      if (i != j) {
	iw[ degree[i] ] = j+1;
	iw[ degree[j] ] = i+1;
	degree[i] ++;
	degree[j] ++;
      }
    }
  }

  taucs_printf("taucs_ccs_amd: calling amd matrix is %dx%d, nnz=%d\n",
	     n,n,nnz);

  if (!strcmp(which,"mmd")) 
    amdexa_(&n, pe, iw, len, &iwlen, &pfree, nv, next,
	    last, head, elen, degree, &ncmpa, w, &iovflo);
  else if (!strcmp(which,"md")) 
    amdtru_(&n, pe, iw, len, &iwlen, &pfree, nv, next,
	  last, head, elen, degree, &ncmpa, w, &iovflo);
  else if (!strcmp(which,"amd")) 
    amdbar_(&n, pe, iw, len, &iwlen, &pfree, nv, next,
	    last, head, elen, degree, &ncmpa, w, &iovflo);
  else {
    taucs_printf("taucs_ccs_amd: WARNING - invalid ordering requested (%s)\n",which);
    return;
  }

  taucs_printf("taucs_ccs_amd: amd returned. optimal iwlen=%d (in this run was %d), %d compressions\n",
	     pfree,iwlen,ncmpa);
  /*
  {
    FILE* f;
    f=fopen("p.ijv","w");
    for (i=0; i<n; i++) fprintf(f,"%d\n",last[i]);
    fclose(f);
  }
  */

  taucs_free(pe    );
  taucs_free(degree);
  taucs_free(nv    );
  taucs_free(next  );
  /* free(last  ); */
  taucs_free(head  );
  taucs_free(elen  );
  taucs_free(w     );
  /* free(len   ); */
  taucs_free(iw    );
  
  for (i=0; i<n; i++) last[i] --;
  for (i=0; i<n; i++) len[ last[i] ] = i;

  *perm    = last;
  *invperm = len;
#endif
}

/*********************************************************/
/* Interface to MMD                                      */
/*********************************************************/

static void
taucs_ccs_genmmd(taucs_ccs_matrix* m, 
		 int** perm, int** invperm,
		 char* which)
{
#ifndef TAUCS_CONFIG_GENMMD
  taucs_printf("taucs_ccs_genmmd: GENMMD routines not linked.\n");
  *perm    = NULL;
  *invperm = NULL;
  return;
#else
  int  n, maxint, delta, nofsub;
  int* xadj;
  int* adjncy;
  int* invp;
  int* prm;
  int* dhead;
  int* qsize;
  int* llist;
  int* marker;

  int* len;
  int* next;

  int  nnz,i,j,ip;
  
  /*taucs_printf("taucs_ccs_genmmd: starting (%s)\n",which);*/

  if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
    taucs_printf("taucs_ccs_genmmd: GENMMD ordering only works on symmetric matrices.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }
  /* this routine may actually work on UPPER as well */
  if (!(m->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_genmmd: the lower part of the matrix must be represented.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }

  *perm    = NULL;
  *invperm = NULL;

  n   = m->n;
  nnz = (m->colptr)[n];
  
  /* I copied the value of delta and the size of */
  /* from SuperLU. Sivan                         */

  delta = 1; /* DELTA is a parameter to allow the choice of nodes
		whose degree <= min-degree + DELTA. */
  delta = 1; /* DELTA is a parameter to allow the choice of nodes
		whose degree <= min-degree + DELTA. */
  /*maxint = 2147483648;*/ /* 2**31-1, for 32-bit only! */
  maxint = 32000;

  assert(sizeof(int) == 4);
  maxint = 2147483647; /* 2**31-1, for 32-bit only! */

  xadj   = (int*) taucs_malloc((n+1)     * sizeof(int));
  adjncy = (int*) taucs_malloc((2*nnz-n) * sizeof(int));
  invp   = (int*) taucs_malloc((n+1)     * sizeof(int));
  prm    = (int*) taucs_malloc(n         * sizeof(int));
  dhead  = (int*) taucs_malloc((n+1)     * sizeof(int));
  qsize  = (int*) taucs_malloc((n+1)     * sizeof(int));
  llist  = (int*) taucs_malloc(n         * sizeof(int));
  marker = (int*) taucs_malloc(n         * sizeof(int));

  if (!xadj || !adjncy || !invp || !prm 
      || !dhead || !qsize || !llist || !marker) {
    taucs_free(xadj  );
    taucs_free(adjncy);
    taucs_free(invp  );
    taucs_free(prm   );
    taucs_free(dhead );
    taucs_free(qsize );
    taucs_free(llist );
    taucs_free(marker);
    return;
  }

  len  = dhead; /* we reuse space */
  next = qsize; /* we reuse space */

  for (i=0; i<n; i++) len[i] = 0;

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      /*i = (m->rowind)[ip] - (m->indshift);*/
      i = (m->rowind)[ip];
      if (i != j) {
	len[i] ++;
	len[j] ++;
      } else {
	/*len[i] ++;*/
      }
    }
  }

  xadj[0] = 1;
  for (i=1; i<=n; i++) xadj[i] = xadj[i-1] + len[i-1];

  /*for (i=0; i<=n; i++) printf("xadj[%d]=%d\n",i,xadj[i]);*/
  
  /* use degree as a temporary */

  for (i=0; i<n; i++) next[i] = xadj[i] - 1;

  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      /*i = (m->rowind)[ip] - (m->indshift);*/
      i = (m->rowind)[ip];
      assert( next[i] < 2*nnz-n );
      assert( next[j] < 2*nnz-n );
      if (i != j) {
	adjncy[ next[i] ] = j+1;
	adjncy[ next[j] ] = i+1;
	next[i] ++;
	next[j] ++;
      } else {
	/*
        adjncy[ next[i] ] = j+1;
	next[i] ++;
	*/
      }
    }
  }

  /*
  for (j=0; j<n; j++) {
    qsort(adjncy + (xadj[j] - 1),
	  xadj[j+1] - xadj[j],
	  sizeof(int),
	  compare_ints);
    printf("+++ %d: ",j+1);
    for (ip=xadj[j]-1; ip<xadj[j+1]-1;ip++)
      printf("%d ",adjncy[ip]);
    printf("\n");
  }
  */

  /*
  taucs_printf("taucs_ccs_genmmd: calling genmmd, matrix is %dx%d, nnz=%d\n",
	     n,n,nnz);
  */

  genmmd_(&n,
	  xadj, adjncy,
	  invp,prm,
	  &delta,
	  dhead,qsize,llist,marker,
	  &maxint,&nofsub);


  /*taucs_printf("taucs_ccs_genmmd: genmmd returned.\n");*/

  /*
  {
    FILE* f;
    f=fopen("p.ijv","w");
    for (i=0; i<n; i++) fprintf(f,"%d %d\n",prm[i],invp[i]);
    fclose(f);
  }
  */

  taucs_free(marker);
  taucs_free(llist );
  taucs_free(qsize );
  taucs_free(dhead );
  taucs_free(xadj  );
  taucs_free(adjncy);
  
  for (i=0; i<n; i++) prm[i] --;
  for (i=0; i<n; i++) invp[ prm[i] ] = i;

  *perm    = prm;
  *invperm = invp;
#endif
}

/*********************************************************/
/* No-fill ordering for trees                            */
/*********************************************************/

static void 
taucs_ccs_treeorder(taucs_ccs_matrix* m,
		    int** perm,
		    int** invperm)
{
  int  n,nnz,i,j,ip,k,p,nleaves;
  int* adjptr;
  int* adj;
  int* len;
  int* ptr;
  int* degree;
  int* leaves;

  if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
    taucs_printf("taucs_ccs_treeorder: tree ordering only works on symmetric matrices.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }
  /* this routine may actually work on UPPER as well */
  if (!(m->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_treeorder: the lower part of the matrix must be represented.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }

  n   = m->n;
  nnz = (m->colptr)[n];
  
  taucs_printf("taucs_ccs_treeorder: starting, matrix is %dx%d, # edges=%d\n",
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
    *perm = *invperm = NULL;
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

  /*
  taucs_printf("taucs_ccs_treeorder: %d initial leaves: ",nleaves);
  for (i=0; i<nleaves; i++) 
    taucs_printf("%d ",leaves[i]);
  taucs_printf("\n");
  */

  for (i=0; i<n; i++) {
    nleaves--;
    if (nleaves <= 0) {
      /* not a tree */
      taucs_free(adj);
      taucs_free(adjptr);
      taucs_free(len);
      taucs_free(leaves);
      taucs_free(degree);
      taucs_free(*perm);
      taucs_free(*invperm);
      *perm = *invperm = NULL;
    }
    j = leaves[nleaves];

    /*taucs_printf("taucs_ccs_treeorder: next leaf is %d, degree=%d\n",j,len[j]);*/
    
    (*perm)   [ i ] = j;
    (*invperm)[ j ] = i;
    
    if (len[j] > 0) {
      if (len[j] != 1) {
	/* not a tree */
	taucs_free(adj);
	taucs_free(adjptr);
	taucs_free(len);
	taucs_free(leaves);
	taucs_free(degree);
	taucs_free(*perm);
	taucs_free(*invperm);
	*perm = *invperm = NULL;
      }
      p = adj[ adjptr[j] ]; 
      
      /*taucs_printf("taucs_ccs_treeorder: parent of %d is %d\n",j,p);*/

      for (k = 0; k < len[p]; k++)
	if (adj[ adjptr[p] + k ] == j) break;

      if ( k >= len[p] ) { /* otherwise j does not show up in p's adjacency list */
	/* not a tree */
	taucs_free(adj);
	taucs_free(adjptr);
	taucs_free(len);
	taucs_free(leaves);
	taucs_free(degree);
	taucs_free(*perm);
	taucs_free(*invperm);
	*perm = *invperm = NULL;
      }

      /* now delete j from p's adjacency list and compress */
      len[p] --;
      for (; k < len[p]; k++)
	adj[ adjptr[p] + k ] = adj[ adjptr[p] + k+1 ];

      if (len[p] == 1) {  /* degree was higher and now is 1 */
	leaves[ nleaves ] = p;
	nleaves++;
      }
    }
  }

  taucs_free(adj);
  taucs_free(adjptr);
  taucs_free(len);
  taucs_free(leaves);
  taucs_free(degree);

  /*
  taucs_printf("taucs_ccs_treeorder: ordering: ");
  for (i=0; i<n; i++) 
    taucs_printf("%d ",(*perm)[i]);
  taucs_printf("\n");
  */

  taucs_printf("taucs_ccs_treeorder: done\n");
}

/*********************************************************/
/* Interface to METIS                                    */
/*********************************************************/

/* from stuct.h in metis */
typedef int idxtype; 
/* from metis.h */
void METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *);

static void 
taucs_ccs_metis(taucs_ccs_matrix* m, 
		int** perm, int** invperm,
		char* which)
{
#ifndef TAUCS_CONFIG_METIS
  taucs_printf("taucs_ccs_metis: METIS routines not linked.\n");
  *perm    = NULL;
  *invperm = NULL;
  return;
#else
  int  n,nnz,i,j,ip;
  int* xadj;
  int* adj;
  int  num_flag     = 0;
  int  options_flag = 0;
  int* len;
  int* ptr;

  /* taucs_printf("taucs_ccs_metis: starting (%s)\n",which); */

  if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
    taucs_printf("taucs_ccs_treeorder: METIS ordering only works on symmetric matrices.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }
  /* this routine may actually work on UPPER as well */
  if (!(m->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_metis: the lower part of the matrix must be represented.\n");
    *perm    = NULL;
    *invperm = NULL;
    return;
  }

  n   = m->n;
  nnz = (m->colptr)[n];
  
  *perm    = (int*) taucs_malloc(n * sizeof(int));
  *invperm = (int*) taucs_malloc(n * sizeof(int));

  xadj = (int*) taucs_malloc((n+1) * sizeof(int));
  /* Change suggested by Yifan Hu for diagonal matrices */
  /* and for matrices with no diagonal */
  /* adj  = (int*) taucs_malloc(2*(nnz-n) * sizeof(int));*/
  adj  = (int*) taucs_malloc(2* nnz * sizeof(int));

  if (!(*perm) || !(*invperm) || !xadj || !adj) {
    taucs_free(*perm);
    taucs_free(*invperm);
    taucs_free(xadj);
    taucs_free(adj);
    *perm = *invperm = NULL;
    return;
  }

  /* assert(*perm && *invperm && xadj && adj);*/

  ptr = len = *perm;

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

  xadj[0] = 0;
  for (i=1; i<=n; i++) xadj[i] = xadj[i-1] + len[i-1];
  
  for (i=0; i<n; i++) ptr[i] = xadj[i];

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

  /* taucs_printf("taucs_ccs_metis: calling metis matrix is %dx%d, nnz=%d\n", */
	     /* n,n,nnz); */

  METIS_NodeND(&n,
	       xadj,adj,
	       &num_flag, &options_flag,
	       *perm,*invperm);

  /* taucs_printf("taucs_ccs_metis: metis returned\n"); */

  /*
  {
    FILE* f;
    f=fopen("p.ijv","w");
    for (i=0; i<n; i++) fprintf(f,"%d\n",last[i]);
    fclose(f);
  }
  */

  taucs_free(xadj);
  taucs_free(adj);
#endif
}

/*********************************************************/
/* RANDOM PERMUTATION                                    */
/*********************************************************/

static void 
taucs_ccs_randomperm(int n,int** perm, int** invperm)
{
  int i;

  *perm    = (int*) taucs_malloc(n * sizeof(int));
  *invperm = (int*) taucs_malloc(n * sizeof(int));
  if (!(*perm) || !(*invperm)) {
    taucs_free(*perm); taucs_free(*invperm);
    *perm = *invperm = NULL;
    taucs_printf("taucs_ccs_randomperm: out of memory for permutation\n");
    return;
  }

  for (i=0; i<n; i++) (*perm)[i] = i;

  for (i=0; i<n; i++) {
    int i1, i2;
    int t;

    i1 = rand() % (n - i);
    i2 = n - i - 1;
    
    t = (*perm)[i1];
    (*perm)[i1] = (*perm)[i2];
    (*perm)[i2] = t;
  }

  for (i=0; i<n; i++) (*invperm)[(*perm)[i]] = i;
  return;
}

/*********************************************************/
/* MAIN ORDERING ROUTINE                                 */
/*********************************************************/

void 
taucs_ccs_order(taucs_ccs_matrix* m, 
		int** perm, int** invperm,
		char* which)
{
  if (!strcmp(which,"mmd") || !strcmp(which,"amd") || !strcmp(which,"md")) 
    taucs_ccs_amd(m,perm,invperm,which);
  else if (!strcmp(which,"metis"))
    taucs_ccs_metis(m,perm,invperm,which);
  else if (!strcmp(which,"genmmd"))
    taucs_ccs_genmmd(m,perm,invperm,which);
  else if (!strcmp(which,"colamd"))
    taucs_ccs_colamd(m,perm,invperm,which);
  else if (!strcmp(which,"random"))
    taucs_ccs_randomperm(m->n,perm,invperm);
  else if (!strcmp(which,"tree")) {
    taucs_ccs_treeorder(m,perm,invperm);
    if (*perm == NULL) /* perhaps the graph of the matrix is not a tree */
      taucs_ccs_metis(m,perm,invperm,"metis");
  }
  else if (!strcmp(which,"identity")) {
    int i;
    *perm    = (int*) taucs_malloc((m->n) * sizeof(int));
    *invperm = (int*) taucs_malloc((m->n) * sizeof(int));
    if (!(*perm) || !(*invperm)) {
      taucs_free(*perm); taucs_free(*invperm);
      *perm = *invperm = NULL;
      taucs_printf("taucs_ccs_order: out of memory for identity permutation\n");
      return;
    }
    for (i=0; i<m->n; i++) (*perm)[i] = (*invperm)[i] = i;
    return;
  }
  else {
    taucs_printf("taucs_ccs_order: invalid ordering requested (%s)\n",which);
    *perm = *invperm = NULL;
  }
}

/*********************************************************/
/*                                                       */
/*********************************************************/

#endif /* TAUCS_CORE_GENERAL */

/*********************************************************/
/*                                                       */
/*********************************************************/

