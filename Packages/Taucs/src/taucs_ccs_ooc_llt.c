
/*
 * taucs_ccs_ooc_llt.c 
 *
 * Out-of-core sparse Cholesky factorization
 *
 * authors: Vladimir Rotkin & Sivan Toledo 
 *
 * Copyright, 2001.
 */

/*************************************************************/
/*                                                           */
/*************************************************************/

#include <stdio.h>
#include <stdlib.h>


/*#include <unistd.h>*/
/*#include <sys/uio.h>*/

#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "taucs.h"

#define FALSE 0
#define TRUE  1

/* #define BLAS_FLOPS_CUTOFF  1000.0 */

#define BLAS_FLOPS_CUTOFF  -1.0
#define SOLVE_DENSE_CUTOFF 5

/* number of matrices in the header of the file */
/*#define IO_BASE    6*/
#define IO_BASE    7
/* multiple files of at most 1 GB or single file */
#define MULTIFILE 0

#ifndef TAUCS_CORE_GENERAL

/*************************************************************/
/* structures                                                */
/*************************************************************/

typedef struct {
  int     sn_size;
  int     n;
  int*    rowind;

  int     up_size;
  int*    sn_vertices;
  int*    up_vertices;
  taucs_datatype* f1;
  taucs_datatype* f2;
  taucs_datatype* u;
} supernodal_frontal_matrix;

typedef struct {
  char    uplo;     /* 'u' for upper, 'l' for lower, ' ' don't know; prefer lower. */
  int     n;        /* size of matrix */
  int     n_sn;     /* number of supernodes */

  int* parent;      /* supernodal elimination tree */
  int* first_child; 
  int* next_child;
  int* ipostorder;
  int* col_to_sn_map;     

  int* sn_size;     /* size of supernodes (diagonal block) */
  int* sn_up_size;  /* size of subdiagonal update blocks   */
  int** sn_struct;  /* row structure of supernodes         */

  taucs_datatype** sn_blocks; /* supernode blocks        */
  taucs_datatype** up_blocks; /* update blocks           */
} supernodal_factor_matrix;


/*************************************************************/
/* for qsort                                                 */
/*************************************************************/
#if 0
static int compare_ints(void* vx, void* vy)
{
  int* ix = (int*)vx;
  int* iy = (int*)vy;
  if (*ix < *iy) return -1;
  if (*ix > *iy) return  1;
  return 0;
}
#endif

static int* compare_indirect_map;
static int compare_indirect_ints(const void* vx, const void* vy)/*(void* vx, void* vy) omer*/
{
  int* ix = (int*)vx;
  int* iy = (int*)vy;
  if (compare_indirect_map[*ix] < compare_indirect_map[*iy]) return -1;
  if (compare_indirect_map[*ix] > compare_indirect_map[*iy]) return  1;
  return 0;
}

/*************************************************************/
/* radix sort                                                */
/*************************************************************/

#if 0

/* NCOUNTS = 2^LOGRADIX */

#define RADIX_SORT_LOGRADIX 4
#define RADIX_SORT_NCOUNTS  16

static unsigned int counts[RADIX_SORT_NCOUNTS];

static int
radix_sort(unsigned int* x, int n)
{
  int i;
  unsigned int mask;

  unsigned int  ncounts;

  unsigned int* y;
  unsigned int* to;
  unsigned int* from;

  unsigned int v;
  unsigned int partialsum;
  unsigned int next;
  unsigned int shift;
  unsigned int bits_sorted;

  if (RADIX_SORT_LOGRADIX >= 8*sizeof(unsigned int)) {
    taucs_printf("radix sort: radix too large.\n");
    /* the computation of ncounts will fail */
    return 0;
  }

  mask    = 0;
  ncounts = 1;
  for (i=0; i<RADIX_SORT_LOGRADIX; i++) {
    mask = (mask << 1) | 1;
    ncounts = ncounts << 1;
  }

  assert(ncounts==RADIX_SORT_NCOUNTS);

  y      = (unsigned int*) taucs_malloc(n       * sizeof(unsigned int));
  if (!y) {
    taucs_printf("radix sort: out of memory.\n");
    return -1;
  }

  from = x;
  to   = y;

  bits_sorted = 0;
  while(bits_sorted < 8*sizeof(unsigned int)) {
    for (i=0; i<ncounts; i++) counts[i] = 0;

    for (i=0; i<n; i++) {
      v = (from[i] >> bits_sorted) & mask;
      assert(v < ncounts);
      counts[v] ++;
    }

    partialsum = 0;
    for (i=0; i<ncounts; i++) {
      /*printf("<%d ",counts[i]);*/
      next = counts[i];
      counts[i] = partialsum;
      /*printf("%d>\n",counts[i]);*/
      partialsum = partialsum + next;
    }

    for (i=0; i<n; i++) {
      v = (from[i] >> bits_sorted) & mask;
      assert(counts[v] < n);
      to[counts[v]] = from[i];
      counts[v] ++;
    }
    /*
    printf("===========\n");
    for (i=0; i<n; i++) printf(">>%d>> %08x\n",bits_sorted,to[i]);
    printf("===========\n");
    */

    bits_sorted += RADIX_SORT_LOGRADIX;
    if (from == x) {
      from = y;
      to   = x;
    } else {
      from = x;
      to   = y;
    } 
  }

  if (from == y) 
    for (i=0; i<n; i++) x[i] = y[i];

  taucs_free(y);

  return 0;
}
#endif

/*************************************************************/
/* create and free the factor object                         */
/*************************************************************/

static supernodal_factor_matrix*
multifrontal_supernodal_create()
{
  supernodal_factor_matrix* L;
  
  L = (supernodal_factor_matrix*) taucs_malloc(sizeof(supernodal_factor_matrix));
  if (!L) return NULL;
  L->uplo      = 'l';
  L->n         = -1; /* unused */

  L->sn_struct   = NULL;
  L->sn_size     = NULL;
  L->sn_up_size  = NULL;
  L->parent      = NULL;
  L->col_to_sn_map = NULL;
  L->first_child = NULL;
  L->next_child  = NULL;
  L->ipostorder  = NULL;
  L->sn_blocks     = NULL;
  L->up_blocks     = NULL;

  return L;
}

static
void ooc_supernodal_factor_free(void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int sn;
  
  taucs_free(L->parent);
  taucs_free(L->first_child);
  taucs_free(L->next_child);
  taucs_free(L->col_to_sn_map);

  taucs_free(L->sn_size);
  taucs_free(L->sn_up_size);
  for (sn=0; sn<L->n_sn; sn++) {
    taucs_free(L->sn_struct[sn]);
    taucs_free(L->sn_blocks[sn]);
    taucs_free(L->up_blocks[sn]);
  }

  taucs_free(L->sn_struct);
  taucs_free(L->sn_blocks);
  taucs_free(L->up_blocks);

  taucs_free(L);
}

static void
recursive_symbolic_elimination(int            j,
			       taucs_ccs_matrix* A,
			       int            first_child[],
			       int            next_child[],
			       int*           n_sn,
			       int            sn_size[],
			       int            sn_up_size[],
			       int*           sn_rowind[],
			       int            sn_first_child[], 
			       int            sn_next_child[], 
			       int            rowind[],
			       int            column_to_sn_map[],
			       int            map[],
			       int            do_order,
			       int            ipostorder[],
			       double         given_mem,
			       void           (*sn_struct_handler)(),
			       void*          sn_struct_handler_arg
			       )
{
  int  i,ip,c,c_sn;
  int  in_previous_sn;
  int  nnz = 0; /* to supress a warning */

  for (c=first_child[j]; c != -1; c = next_child[c]) {
    recursive_symbolic_elimination(c,A,
				   first_child,next_child,
				   n_sn,
				   sn_size,sn_up_size,sn_rowind,
				   sn_first_child,sn_next_child,
				   rowind, /* temporary */
				   column_to_sn_map,
				   map,
				   do_order,ipostorder,given_mem,
				   sn_struct_handler,sn_struct_handler_arg
				   );
  }

  
  in_previous_sn = 1;
  if (j == A->n) 
    in_previous_sn = 0; /* this is not a real column */
  else if (first_child[j] == -1) 
    in_previous_sn = 0; /* this is a leaf */
  else if (next_child[first_child[j]] != -1) 
    in_previous_sn = 0; /* more than 1 child */
  else if ((double)sn_up_size[column_to_sn_map[first_child[j]]]
	   *(double)(sn_size[column_to_sn_map[first_child[j]]]+1)
	   *sizeof(taucs_datatype) > given_mem)
    in_previous_sn = 0; /* size of supernode great than given memory */
  else { 
    /* check that the structure is nested */
    /* map contains child markers         */

    c=first_child[j];
    for (ip=(A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      in_previous_sn = in_previous_sn && (map[i] == c);
    }
  }

  if (in_previous_sn) {
    c = first_child[j];
    c_sn = column_to_sn_map[c];
    column_to_sn_map[j] = c_sn;

    /* swap row indices so j is at the end of the */
    /* supernode, not in the update indices       */
    for (ip=sn_size[c_sn]; ip<sn_up_size[c_sn]; ip++) 
      if (sn_rowind[c_sn][ip] == j) break;
    assert(ip<sn_up_size[c_sn]);
    sn_rowind[c_sn][ip] = sn_rowind[c_sn][sn_size[c_sn]];
    sn_rowind[c_sn][sn_size[c_sn]] = j;

    /* mark the nonzeros in the map */
    for (ip=sn_size[c_sn]; ip<sn_up_size[c_sn]; ip++) 
      map[ sn_rowind[c_sn][ip] ] = j;

    sn_size   [c_sn]++;
    if((double)sn_size[c_sn]*(double)sn_up_size[c_sn]>=(1024.0*1024.0*1024.0))
      taucs_printf("debug!!!: sn_size[%d] = %d sn_up_size[%d] = %d\n ",c_sn,sn_size[c_sn],c_sn,sn_up_size[c_sn]);
    /* return c_sn; */
    return;
  }

  /* we are in a new supernode */
  if (j < A->n) {
    nnz = 1;
    rowind[0] = j;
    map[j]    = j;
    
    for (c=first_child[j]; c != -1; c = next_child[c]) {
      c_sn = column_to_sn_map[c];
      for (ip=sn_size[c_sn]; ip<sn_up_size[c_sn]; ip++) {
	i = sn_rowind[c_sn][ip];
	if (i > j && map[i] != j) { /* new row index */
	  map[i] = j;
	  rowind[nnz] = i;
	  nnz++;
	}
      }
      if (sn_struct_handler)
	(*sn_struct_handler)(sn_struct_handler_arg,
			     c_sn,sn_up_size[c_sn],&(sn_rowind[c_sn]));
    }

    for (ip=(A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (map[i] != j) { /* new row index */
	map[i] = j;
	rowind[nnz] = i;
	nnz++;
      }
    }
  }

  /* append childs from root*/
  if (j == A->n) {
    for (c=first_child[j]; c != -1; c = next_child[c]) {
      c_sn = column_to_sn_map[c];
      if (sn_struct_handler)
	(*sn_struct_handler)(sn_struct_handler_arg,
			     c_sn,sn_up_size[c_sn],&(sn_rowind[c_sn]));
    }
  }

  /*printf("children of sn %d: ",*n_sn);*/
  for (c=first_child[j]; c != -1; c = next_child[c]) {
    c_sn = column_to_sn_map[c];
    /*printf("%d ",c_sn);*/
    if (c==first_child[j])
      sn_first_child[*n_sn] = c_sn;
    else {
      sn_next_child[ c_sn ] = sn_first_child[*n_sn];
      sn_first_child[*n_sn] = c_sn;
    }
  }
  /*printf("\n");*/

  if (j < A->n) {
    column_to_sn_map[j] = *n_sn;
    sn_size   [*n_sn] = 1;
    sn_up_size[*n_sn] = nnz;
    sn_rowind [*n_sn] = (int*) taucs_malloc(nnz * sizeof(int));
    for (ip=0; ip<nnz; ip++) sn_rowind[*n_sn][ip] = rowind[ip];
    if (do_order) {
      /* Sivan and Vladimir: we think that we can sort in */
      /* column order, not only in etree postorder.       */
      /*
	radix_sort(sn_rowind [*n_sn],nnz);
	qsort(sn_rowind [*n_sn],nnz,sizeof(int),compare_ints);
      */
      compare_indirect_map = ipostorder;
      qsort(sn_rowind [*n_sn],nnz,sizeof(int),compare_indirect_ints);
    }
    assert(sn_rowind [*n_sn][0] == j);
    (*n_sn)++;
  }
}

static
void recursive_postorder(int  j,
			 int  first_child[],
			 int  next_child[],
			 int  postorder[],
			 int  ipostorder[],
			 int* next
			 )
{
  int c;


  for (c=first_child[j]; c != -1; c = next_child[c]) {
    /*printf("*** %d is child of %d\n",c,j);*/
    recursive_postorder(c,first_child,next_child,
			postorder,ipostorder,next
			);
  }
  /*  printf(">>> j=%d next=%d\n",j,*next);*/
  if (postorder)  postorder [*next] = j;
  if (ipostorder) ipostorder[j] = *next;
  (*next)++;
}

static void
taucs_ccs_ooc_symbolic_elimination(taucs_ccs_matrix* A,
				   void* vL,
				   int do_order,
				   int do_column_to_sn_map,
				   double given_mem,
				   void           (*sn_struct_handler)(),
				   void*          sn_struct_handler_arg
				   )
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int* first_child;
  int* next_child;
  int j;
  int* column_to_sn_map;
  int* map;
  int* rowind;
  int* parent;
  int* ipostorder;

  L->n         = A->n;
  L->sn_struct = (int**)taucs_malloc((A->n  )*sizeof(int*));
  L->sn_size   = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->sn_up_size   = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->first_child = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->next_child  = (int*) taucs_malloc((A->n+1)*sizeof(int));

  column_to_sn_map = (int*)taucs_malloc((A->n+1)*sizeof(int));
  map              = (int*) taucs_malloc((A->n+1)*sizeof(int));

  first_child = (int*) taucs_malloc(((A->n)+1)*sizeof(int));
  next_child  = (int*) taucs_malloc(((A->n)+1)*sizeof(int));
    
  rowind      = (int*) taucs_malloc((A->n)*sizeof(int));

  taucs_printf("STARTING SYMB 1\n");

  /* compute the vertex elimination tree */
  parent      = (int*)taucs_malloc((A->n+1)*sizeof(int));
  taucs_ccs_etree(A,parent,NULL,NULL,NULL);
  for (j=0; j <= (A->n); j++) first_child[j] = -1;
  for (j = (A->n)-1; j >= 0; j--) {
    int p = parent[j];
    next_child[j] = first_child[p];
    first_child[p] = j;
  }
  taucs_free(parent);

  taucs_printf("STARTING SYMB 2\n");

  ipostorder = (int*)taucs_malloc((A->n+1)*sizeof(int));
  { 
    int next = 0;
    /*int* postorder = (int*)taucs_malloc((A->n+1)*sizeof(int));*/
    recursive_postorder(A->n,first_child,next_child,
			NULL,
			ipostorder,&next);
    /*
    printf("ipostorder ");
    for (j=0; j <= (A->n); j++) printf("%d ",ipostorder[j]);
    printf("\n");
    printf(" postorder ");
    for (j=0; j <= (A->n); j++) printf("%d ",postorder[j]);
    printf("\n");
    */
  }

  taucs_printf("STARTING SYMB 3\n");

  L->n_sn = 0;
  for (j=0; j < (A->n); j++) map[j] = -1;
  for (j=0; j <= (A->n); j++) (L->first_child)[j] = (L->next_child)[j] = -1;

  taucs_printf("STARTING SYMB\n");

  recursive_symbolic_elimination(A->n,
				 A,
				 first_child,next_child,
				 &(L->n_sn),
				 L->sn_size,L->sn_up_size,L->sn_struct,
				 L->first_child,L->next_child,
				 rowind,
				 column_to_sn_map,
				 map,
				 do_order,ipostorder,given_mem,
				 sn_struct_handler,sn_struct_handler_arg
				 );

  taucs_printf("AFTER SYMB\n");

  {
    double nnz   = 0.0;
    double flops = 0.0;
    int sn,i,colnnz;
    for (sn=0; sn<(L->n_sn); sn++) {
      for (i=0, colnnz = (L->sn_up_size)[sn]; 
	   i<(L->sn_size)[sn]; 
	   i++, colnnz--) {
	flops += 1.0 + ((double)(colnnz)) * ((double)(colnnz));
	nnz   += (double) (colnnz);
      }
    }
    taucs_printf("\t\tSymbolic Analysis of LL^T: %.2e nonzeros, %.2e flops\n",
		 nnz, flops);
  }

  /*
  {
    int i;
    printf("c2sn: ");
    for (i=0; i<A->n; i++) printf("%d ",column_to_sn_map[i]);
    printf("\n");
  }
  */
  
  L->sn_struct = (int**)taucs_realloc( L->sn_struct,(L->n_sn  )*sizeof(int*));
  L->sn_size   = (int*) taucs_realloc( L->sn_size,(L->n_sn+1)*sizeof(int));
  L->sn_up_size   = (int*) taucs_realloc(L->sn_up_size,(L->n_sn+1)*sizeof(int));
  L->first_child = (int*) taucs_realloc(L->first_child,(L->n_sn+1)*sizeof(int));
  L->next_child  = (int*) taucs_realloc(L->next_child,(L->n_sn+1)*sizeof(int));

  L->sn_blocks     = taucs_calloc((L->n_sn), sizeof(taucs_datatype*)); /* so we can free before allocation */
  L->up_blocks     = taucs_calloc((L->n_sn), sizeof(taucs_datatype*));

  taucs_free(rowind);
  taucs_free(map);

  if(do_column_to_sn_map)
    L->col_to_sn_map = column_to_sn_map;
  else
    taucs_free(column_to_sn_map);

  taucs_free(next_child);
  taucs_free(first_child);
  taucs_free(ipostorder);
}

/*************************************************************/
/* left-looking factor routines                              */
/*************************************************************/

static void
recursive_leftlooking_supernodal_update(int J,int K,
					int bitmap[],
					taucs_datatype* dense_update_matrix,
					taucs_ccs_matrix* A,
					supernodal_factor_matrix* L)
{
  int i,j,ir;
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  int sn_size_father = (L->sn_size)[J];
  int sn_up_size_father = (L->sn_up_size)[J];
  int sn_size_child = (L->sn_size)[K];
  int sn_up_size_child = (L->sn_up_size)[K];
  int exist_upd=0;
  int first_row = 0; /* to suppress a warning */
  int row_count=0;
  int PK,M,N,LDA,LDB,LDC;

  for(i=0;i<sn_size_father;i++) {
    bitmap[L->sn_struct[J][i]]=i+1;
  }

  for(i=sn_size_father;i<sn_up_size_father;i++)
    bitmap[L->sn_struct[J][i]] = i - sn_size_father + 1;

  for(i=sn_size_child;i<sn_up_size_child;i++)
    /* is this row index included in the columns of sn J? */
    if(bitmap[L->sn_struct[K][i]]
       && L->sn_struct[K][i] <= L->sn_struct[J][sn_size_father-1]) {
      if(!exist_upd) first_row = i;
      row_count++;
      exist_upd = 1;
      /*taucs_printf("update from K = %d to J = %d \n",K,J);*/
      /* loop over columns of sn K */
            
      /* for(j=0;j<sn_size_child;j++)
	for(ir=i;ir<sn_up_size_child;ir++)
	  if( L->sn_struct[K][ir] <= L->sn_struct[J][sn_size_father-1]){
	    L->sn_blocks[J][ (bitmap[L->sn_struct[K][i]]-1)*(L->sn_blocks_ld[J])+(bitmap[L->sn_struct[K][ir]]-1)] -= L->up_blocks[K][j*(L->up_blocks_ld[K])+ir-sn_size_child]* L->up_blocks[K][j*L->up_blocks_ld[K]+i-sn_size_child];
	    taucs_printf("sn_block: L[%d,%d] = %lf\n",(bitmap[L->sn_struct[K][ir]]-1),(bitmap[L->sn_struct[K][i]]-1),L->sn_blocks[J][ (bitmap[L->sn_struct[K][i]]-1)*(L->sn_blocks_ld[J])+(bitmap[L->sn_struct[K][ir]]-1)]);}
	  else{
	    L->up_blocks[J][ (bitmap[L->sn_struct[K][i]]-1)*(L->up_blocks_ld[J])+(bitmap[L->sn_struct[K][ir]]-1)] -=  L->up_blocks[K][j*L->up_blocks_ld[K]+ir-sn_size_child]* L->up_blocks[K][j*L->up_blocks_ld[K]+i-sn_size_child];
	   taucs_printf("up_block: L[%d,%d] = %lf\n",(bitmap[L->sn_struct[K][ir]]-1),(bitmap[L->sn_struct[K][i]]-1),L->up_blocks[J][ (bitmap[L->sn_struct[K][i]]-1)*(L->up_blocks_ld[J])+(bitmap[L->sn_struct[K][ir]]-1)]);
	   }*/
        }

  if(exist_upd){
    LDA = LDB = (L->sn_up_size)[K]-(L->sn_size)[K];
    M  = sn_up_size_child - first_row ; /* +-1 ? */    
    LDC =  sn_up_size_father;
    N  = row_count; 
    PK = L->sn_size[K];    
   
    /* The GEMM code computes on the upper triangle of the trapezoidal
       matrix, which is junk. */
    /*
    taucs_gemm ("No Conjugate",
		"Conjugate",
		&M,&N,&PK,
		&taucs_one_const,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		&taucs_zero_const,
		dense_update_matrix,&LDC);
    */

    /* This is the HERK+GEMM fix by Elad */
    taucs_herk ("Lower",
		"No Conjugate",
		&N,&PK,
		&taucs_one_real_const,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		&taucs_zero_real_const,
		dense_update_matrix,&LDC);

    if(M-N > 0)
    {
        int newM = M - N;
   
        taucs_gemm ("No Conjugate",
		"Conjugate",
		&newM,&N,&PK,
		&taucs_one_const,
		&(L->up_blocks[K][first_row-sn_size_child+N]),&LDA,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		&taucs_zero_const,
		dense_update_matrix+N,&LDC);
    }
    /* end of GEMM/HERK+GEMM fix */ 

    /*for(j=0;j<row_count;j++)
       for(ir=0;ir<sn_up_size_father;ir++)
	 taucs_printf("dense[%d,%d] = %lf\n",ir,j,dense_update_matrix[j*LDC+ir]);
    */
    for(j=0;j<row_count;j++)
      for(ir=j;ir<row_count;ir++){
	/*
	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] -= dense_update_matrix[j*LDC+ir];
	*/
	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] =
	  taucs_sub(L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] , dense_update_matrix[j*LDC+ir]);

      }

    for(j=0;j<row_count;j++)
      for(ir=row_count;ir<M;ir++){
	/*
	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[J]-L->sn_size[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] -= dense_update_matrix[j*LDC+ir];
	*/
	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[J]-L->sn_size[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] =
	  taucs_sub(L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[J]-L->sn_size[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] , dense_update_matrix[j*LDC+ir]);
	}
    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;
    
    for (child = first_child[K]; child != -1; child = next_child[child]) {
      recursive_leftlooking_supernodal_update(J,child,
					      bitmap,dense_update_matrix,
					      A,L);
    }
  }
  else
    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;

}


static int
leftlooking_supernodal_front_factor(int sn,
				    int* indmap,
				    taucs_ccs_matrix* A,
				    supernodal_factor_matrix* L)
{
  int ip,jp;
  int*    ind;
  taucs_datatype* re;
  int INFO;
  int sn_size = (L->sn_size)[sn];
  int up_size = (L->sn_up_size)[sn] - (L->sn_size)[sn];

  /* creating transform for real indices */
  for(ip=0;ip<(L->sn_up_size)[sn];ip++) indmap[(L->sn_struct)[sn][ip]] = ip;

  for(jp=0;jp<sn_size;jp++) {
    ind = &(A->rowind[A->colptr[ (L->sn_struct)[sn][jp] ]]);
    re  = &(A->taucs_values[A->colptr[ (L->sn_struct)[sn][jp] ]]); 
    for(ip=0;
	ip < A->colptr[ (L->sn_struct)[sn][jp] + 1 ] 
           - A->colptr[ (L->sn_struct)[sn][jp] ];
	ip++) {
      if (indmap[ind[ip]] < sn_size)
	(L->sn_blocks)[sn][sn_size*jp + indmap[ind[ip]]] =
	  taucs_add( (L->sn_blocks)[sn][sn_size*jp + indmap[ind[ip]]] , re[ip]);
      else
	(L->up_blocks)[sn][up_size*jp + indmap[ind[ip]] - sn_size] =
	  taucs_add((L->up_blocks)[sn][up_size*jp + indmap[ind[ip]] - sn_size] , re[ip]);
    }
  }

  /* we use the BLAS through the Fortran interface */

  /* solving of lower triangular system for L */
  if (sn_size)
    taucs_potrf ("LOWER",
		 &sn_size,
		 (L->sn_blocks)[sn],&sn_size,
		 &INFO);

  if (INFO) {
    taucs_printf("\t\tLL^T Factorization: Matrix is not positive definite.\n");
    taucs_printf("\t\t in sn = %d   nonpositive pivot in column %d\n",
		 sn,(L->sn_struct)[sn][INFO-1]);
    return -1;
  }

  /* getting completion for found columns of L */
  if (up_size && sn_size)
    taucs_trsm ("Right",
		"Lower",
		"Conjugate",
		"No unit diagonal",
		&up_size,&sn_size,
		&taucs_one_const,
		(L->sn_blocks)[sn],&sn_size,
		(L->up_blocks)[sn],&up_size);
  
  /* zeroes map */
  for(ip=0;ip<(L->sn_up_size)[sn];ip++) indmap[(L->sn_struct)[sn][ip]] = 0;

  return 0;
}

static int
recursive_leftlooking_supernodal_factor_llt(int sn,       /* this supernode */
					    int is_root,  /* is v the root? */
					    int* map,
					    taucs_ccs_matrix* A,
					    supernodal_factor_matrix* L)
{
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  taucs_datatype* dense_update_matrix = NULL;

  /*  if (!is_root)
      sn_size = L->sn_size[sn];
      else
      sn_size = -1;

      if (!is_root) { 
      (L->sn_blocks   )[sn] = (double*)taucs_calloc(((L->sn_size)[sn])*((L->sn_size)[sn]),sizeof(double));
      
      (L->up_blocks   )[sn] = (double*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn])*((L->sn_size)[sn]),sizeof(double));

    }*/

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if (recursive_leftlooking_supernodal_factor_llt(child,
						    FALSE,
						    map,
						    A,L)) {
      /* failure */
      return -1;
    }
  }    

  if (!is_root) {
    (L->sn_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_size)[sn])*((L->sn_size)[sn]),
						    sizeof(taucs_datatype));
    (L->up_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn])*((L->sn_size)[sn]),
					    sizeof(taucs_datatype));
    if (!dense_update_matrix) 
      dense_update_matrix = (taucs_datatype*) taucs_calloc((L->sn_up_size)[sn]*(L->sn_size)[sn],
						     sizeof(taucs_datatype));
    for (child = first_child[sn]; child != -1; child = next_child[child])
      recursive_leftlooking_supernodal_update(sn,child,
					      map,dense_update_matrix,
					      A,L);
    taucs_free(dense_update_matrix);
    if (leftlooking_supernodal_front_factor(sn,
					    map,
					    A,
					    L)) {
      /* nonpositive pivot */
      return -1;
    }
  }

  return 0;
}

#if 0
void* taucs_ccs_factor_llt_ll(taucs_ccs_matrix* A)
{
  supernodal_factor_matrix* L;
  int i,j,ip,jp;
  int sn,p;
  int* map;
  double wtime, ctime;
 
  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();

  taucs_ccs_ooc_symbolic_elimination(A,L,
				     TRUE /* sort row indices */,
				     FALSE /* don't return col_tosn_map */,
				     1.0/0.0,
				     NULL,NULL /* sn_struct handler*/ );

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);


  map  = (int*)taucs_malloc((A->n+1)*sizeof(int));

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  if (recursive_leftlooking_supernodal_factor_llt((L->n_sn),  
						  TRUE, 
						  map,
						  A,L)) {
    ooc_supernodal_factor_free(L);
    taucs_free(map);

    return NULL;
  }

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSupernodal Left-Looking LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  taucs_free(map);

  return (void*) L;
}
#endif

/*************************************************************/
/* left-looking ooc factor routines                          */
/*************************************************************/

static
void
ooc_sn_struct_handler(void* argument, 
		      int sn, int sn_up_size, int* sn_struct_ptr[])
{
  taucs_io_handle* handle = (taucs_io_handle*) argument;
  taucs_io_append(handle,
		  IO_BASE+sn, /* matrix written in postorder */
		  1,
		  sn_up_size,
		  TAUCS_INT,
		  *sn_struct_ptr);
  taucs_free(*sn_struct_ptr);
  *sn_struct_ptr = NULL;
}

static void
recursive_leftlooking_supernodal_update_ooc(int J,int K,
					    int bitmap[],
					    taucs_datatype* dense_update_matrix,
					    taucs_io_handle* handle,
					    taucs_ccs_matrix* A,
					    supernodal_factor_matrix* L)
{
  int i,j,ir;
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  int sn_size_father = (L->sn_size)[J];
  int sn_up_size_father = (L->sn_up_size)[J];
  int sn_size_child = (L->sn_size)[K];
  int sn_up_size_child = (L->sn_up_size)[K];
  int exist_upd=0;
  int first_row = 0; /* to supress a warning */
  int row_count=0;
  int PK,M,N,LDA,LDB,LDC;

  if(L->sn_up_size[K]-L->sn_size[K]>0){
  
    for(i=0;i<sn_size_father;i++) {
      bitmap[L->sn_struct[J][i]]=i+1;
    }

    for(i=sn_size_father;i<sn_up_size_father;i++){
      bitmap[L->sn_struct[J][i]] = i - sn_size_father + 1;
    }

    L->sn_struct[K] = (int*)taucs_malloc(sn_up_size_child*sizeof(int));
    taucs_io_read(handle,IO_BASE+K,1,sn_up_size_child,TAUCS_INT,L->sn_struct[K]);
    
    for(i=sn_size_child;i<sn_up_size_child;i++){
      /* is this row index included in the columns of sn J? */
      if(bitmap[L->sn_struct[K][i]]
	 && L->sn_struct[K][i] <= L->sn_struct[J][sn_size_father-1]) {
	if(!exist_upd) first_row = i;
	row_count++;
	exist_upd = 1;
      }
    }
  }

  if(exist_upd){
    (L->up_blocks   )[K] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[K]-(L->sn_size)[K])
						   *((L->sn_size)[K]),
						   sizeof(taucs_datatype));
    taucs_io_read(handle,IO_BASE+L->n_sn+2*K+1,
		  (L->sn_up_size)[K]-(L->sn_size)[K],
		  (L->sn_size)[K] ,
		  TAUCS_CORE_DATATYPE,(L->up_blocks)[K]);
   
    LDA = LDB = L->sn_up_size[K]-L->sn_size[K];
    M  = sn_up_size_child - first_row ; /* +-1 ? */    
    LDC =  sn_up_size_father;
    N  = row_count; 
    PK = L->sn_size[K];    
   
    /* The GEMM code computes on the upper triangle of the trapezoidal
       matrix, which is junk. */
    /*
    taucs_gemm ("No Conjugate",
		"Conjugate",
		&M,&N,&PK,
		&taucs_one_const,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		&taucs_zero_const,
		dense_update_matrix,&LDC);
    */

    /* This is the HERK+GEMM fix by Elad */
    taucs_herk ("Lower",
		"No Conjugate",
		&N,&PK,
		&taucs_one_real_const,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		&taucs_zero_real_const,
		dense_update_matrix,&LDC);

    if(M-N > 0)
    {
        int newM = M - N;
   
        taucs_gemm ("No Conjugate",
		"Conjugate",
		&newM,&N,&PK,
		&taucs_one_const,
		&(L->up_blocks[K][first_row-sn_size_child+N]),&LDA,
		&(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		&taucs_zero_const,
		dense_update_matrix+N,&LDC);
    }
    /* end of GEMM/HERK+GEMM fix */ 

    taucs_free((L->up_blocks   )[K]);
    (L->up_blocks   )[K] = NULL;

    for(j=0;j<row_count;j++)
      for(ir=j;ir<row_count;ir++){
	/*
	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] -= dense_update_matrix[j*LDC+ir];
	*/
	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] =
	  taucs_sub(L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] , dense_update_matrix[j*LDC+ir]);
      }

    for(j=0;j<row_count;j++)
      for(ir=row_count;ir<M;ir++){
	/*
	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*((L->sn_up_size)[J]-(L->sn_size)[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] -= dense_update_matrix[j*LDC+ir];
	*/
	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*((L->sn_up_size)[J]-(L->sn_size)[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] =
	  taucs_sub(L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*((L->sn_up_size)[J]-(L->sn_size)[J])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] , dense_update_matrix[j*LDC+ir]);
      }

    taucs_free( L->sn_struct[K]);
    L->sn_struct[K] = NULL;

    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;
    
    for (child = first_child[K]; child != -1; child = next_child[child]) {
      recursive_leftlooking_supernodal_update_ooc(J,child,
						  bitmap,
						  dense_update_matrix,
						  handle,A,L);
    }
  } else {
    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;

    taucs_free( L->sn_struct[K]);
    L->sn_struct[K] = NULL;
  }
}

static int
recursive_append_L(int sn,   /* this supernode */
		   int is_root,/* is v the root?*/
		   taucs_io_handle* handle,
		   supernodal_factor_matrix* L)
{
  int  child;
  /*int  sn_size;*/
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if (recursive_append_L(child,FALSE,handle,L)) {
      /* failure */
      return -1;
    }
  }
    
  if (!is_root) { 
    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn,
		    L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->sn_blocks[sn]);

    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn+1,
		    L->sn_up_size[sn] - L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->up_blocks[sn]); 
    taucs_free((L->sn_blocks   )[sn]);
    taucs_free((L->up_blocks   )[sn]);
    taucs_free((L->sn_struct   )[sn]);
    (L->sn_blocks   )[sn] = NULL;
    (L->up_blocks   )[sn] = NULL;
    (L->sn_struct   )[sn] = NULL;
 }

  return 0;
}

static int
recursive_read_L(int sn,   /* this supernode */
		 int is_root,/* is v the root?*/
		 taucs_io_handle* handle,
		 supernodal_factor_matrix* L)
{
  int  child;
  /*int  sn_size;*/
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if (recursive_read_L(child,FALSE,handle,L)) {
      /* failure */
      return -1;
    }
  }
    
  if (!is_root) { 
    (L->sn_blocks)[sn] = (taucs_datatype*)taucs_calloc(((L->sn_size)[sn])
						 *((L->sn_size)[sn]),
						 sizeof(taucs_datatype));
    (L->up_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn])
						    *((L->sn_size)[sn]),
						    sizeof(taucs_datatype));

    taucs_io_read(handle,IO_BASE+L->n_sn+2*sn,
		  L->sn_size[sn],L->sn_size[sn],
		  TAUCS_CORE_DATATYPE,L->sn_blocks[sn]);

    taucs_io_read(handle,IO_BASE+L->n_sn+2*sn+1,
		  L->sn_up_size[sn] - L->sn_size[sn],L->sn_size[sn],
		  TAUCS_CORE_DATATYPE,L->up_blocks[sn]); 
 }

  return 0;
}

static int
recursive_read_L_cols(int sn,   /* this supernode */
		      int is_root,/* is v the root?*/
		      taucs_io_handle* handle,
		      supernodal_factor_matrix* L)
{
  int  child;
  /*int  sn_size;*/
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if (recursive_read_L_cols(child,FALSE,handle,L)) {
      /* failure */
      return -1;
    }
  }
    
  if (!is_root) { 
    L->sn_struct[sn] = (int*)taucs_malloc((L->sn_up_size)[sn]*sizeof(int));
    taucs_io_read(handle,IO_BASE+sn,1,(L->sn_up_size)[sn],TAUCS_INT,L->sn_struct[sn]);
      
  }

  return 0;
}

static int
recursive_leftlooking_supernodal_factor_llt_ooc(int sn,    /* this supernode */
						int is_root,/* is v the root?*/
						int* map,
						int* sn_in_core,
						taucs_io_handle* handle,
						taucs_ccs_matrix* A,
						supernodal_factor_matrix* L)
{
  int  child;
  /*int  sn_size;*/
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  taucs_datatype* dense_update_matrix = NULL;
  
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if(sn_in_core[child]){

      if (recursive_read_L_cols(child,
				FALSE,
				handle,
				L)) {
	return -1;
      }
      if (recursive_leftlooking_supernodal_factor_llt(child,
						      FALSE,
						      map,
						      A,L)) {
	/* failure */
	return -1;
      }
      if (recursive_append_L(child,
			     FALSE,
			     handle,
			     L)) {
	return -1;
      }

    }
    else
      if (recursive_leftlooking_supernodal_factor_llt_ooc(child,
							  FALSE,
							  map,
							  sn_in_core,
							  handle,
							  A,L)) {
	/* failure */
	return -1;
      }
  }

  if (!is_root) { 
    (L->sn_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_size)[sn])
						    *((L->sn_size)[sn]),
						    sizeof(taucs_datatype));
    
    if(L->sn_up_size[sn]-L->sn_size[sn]>0) (L->up_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn]) *((L->sn_size)[sn]),sizeof(taucs_datatype));

    L->sn_struct[sn] = (int*)taucs_malloc((L->sn_up_size)[sn]*sizeof(int));
    
    taucs_io_read(handle,IO_BASE+sn,1,(L->sn_up_size)[sn],TAUCS_INT,L->sn_struct[sn]);
    
    if (!dense_update_matrix) 
      dense_update_matrix = (taucs_datatype*) taucs_calloc((L->sn_up_size)[sn]*(L->sn_size)[sn],sizeof(taucs_datatype));
    for (child = first_child[sn]; child != -1; child = next_child[child]) {
      recursive_leftlooking_supernodal_update_ooc(sn,child,
						  map,
						  dense_update_matrix,
						  handle,A,L);
    }

    taucs_free(dense_update_matrix);
    if (leftlooking_supernodal_front_factor(sn,
					    map,
					    A,
					    L)) {
      /* nonpositive pivot */
      return -1;
    }
    
    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn,
		    L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->sn_blocks[sn]);
    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn+1,
		    L->sn_up_size[sn] - L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->up_blocks[sn]);
    
    taucs_free((L->sn_blocks)[sn]);
    taucs_free((L->up_blocks)[sn]); 
    taucs_free((L->sn_struct)[sn]);
    (L->sn_blocks   )[sn] = NULL;
    (L->up_blocks   )[sn] = NULL;
    (L->sn_struct   )[sn] = NULL;
  }

  return 0;
}

static double
recursive_compute_supernodes_ll_in_core(int sn,       /* this supernode */
					int is_root,  /* is v the root? */
					double avail_mem,
					int* sn_in_core,
					supernodal_factor_matrix* L)
{
  int  child;
  /*double curr_avail_mem = avail_mem;*/
  double child_mem = 0.0;
  double children_mem = 0.0;
  double total_mem = 0.0;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  
  /*new_panel = 0;*/
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    child_mem=recursive_compute_supernodes_ll_in_core(child,
						      FALSE,
						      avail_mem,
						      sn_in_core,
						      L);
    /*if (child_mem == 0) new_panel = 1; */
    children_mem += child_mem;
  }

  /*if (first_child[sn] == -1) 
    total_mem = (double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype)+(double)(L->sn_up_size)[sn]*sizeof(int);
  else
    total_mem = children_mem + 2*(double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype)+(double)(L->sn_up_size)[sn]*sizeof(int);
  
  if (total_mem <= avail_mem || first_child[sn] == -1) 
    sn_in_core[sn] = 1;
   else 
     sn_in_core[sn] = 0;
  return total_mem;
  */
  total_mem = children_mem + (double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype)+(double)(L->sn_up_size)[sn]*sizeof(int);

  if ((total_mem + (double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype)) <= avail_mem || first_child[sn] == -1) {
    sn_in_core[sn] = 1;
    return total_mem;
  } else {
    sn_in_core[sn] = 0;
    return (total_mem + (double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype));
    }
  
}

static int
recursive_compute_supernodes_in_core_old(int sn,       /* this supernode */
				     int is_root,  /* is v the root? */
				     int avail_mem,
				     int* sn_in_core,
				     supernodal_factor_matrix* L)
{
  int  child;
  int curr_avail_mem = avail_mem;
  int catch_mem = 0;
  int child_mem = 0;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  /*taucs_datatype* dense_update_matrix = NULL;*/


  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    child_mem=recursive_compute_supernodes_in_core_old(child,
						   FALSE,
						   curr_avail_mem,
						   sn_in_core,
						   L);
    if ((int)(child_mem+(L->sn_size)[child]*(L->sn_up_size)[child]*sizeof(taucs_datatype))<curr_avail_mem) sn_in_core[child] = 1;

    catch_mem += child_mem;
  }
  if (!is_root) catch_mem += (L->sn_size)[sn]*(L->sn_up_size)[sn]*sizeof(taucs_datatype);
    return catch_mem;
}

static void
recursive_compute_supernodes_ipostorder(int sn,       /* this supernode */
					int is_root,  /* is v the root? */
					int* current_index_ptr,
					supernodal_factor_matrix* L)
{
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
 
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    recursive_compute_supernodes_ipostorder(child,
					     FALSE,
					     current_index_ptr,
					     L);
  
  }
  L->ipostorder[sn] = *current_index_ptr;
  (*current_index_ptr)++;
 }

/* no panelization */
#if 0
void* taucs_ccs_factor_llt_ll_ooc(taucs_ccs_matrix* A, 
				  double memory)
{
  supernodal_factor_matrix* L;
  int i,j,ip,jp;
  int sn,p;
  int* map;
  int* sn_in_core;
  double wtime, ctime;
  taucs_io_handle* handle;
  int current_index=0;
  double memory_overhead;

  /* compute fixed memory overhead */
  memory_overhead = 
    4.0*(double)((A->n)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((A->n)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((A->n)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((A->n)*sizeof(int));  /* singlefile matrix arrays */  
 
  taucs_printf("max memory overhead %.0f memory %.0f\n",memory_overhead,memory);

  if ( memory - memory_overhead < 
       2.0*(double)((A->n)*sizeof(taucs_datatype)) + 
       2.0*(double)((A->n)*sizeof(int)) ) {
    taucs_printf("\t\ttaucs_ccs_factor_llt_ll_ooc: not enough memory\n");
    return NULL;
  }

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();
  if(!MULTIFILE)
    handle = taucs_io_create_singlefile("/tmp/taucs.L");
  else
    handle = taucs_io_create_multifile("/tmp/taucs.L");

  taucs_io_append(handle,5,1,1,TAUCS_INT,&(A->n));

  taucs_ccs_ooc_symbolic_elimination(A,L,
				     TRUE /* sort row indices */,
				     FALSE /* don't return col_to_sn_map */,
				     (memory - memory_overhead)/3.0,
				     ooc_sn_struct_handler,handle);
  
  /* we now compute an exact memory overhead bound using n_sn */
  memory_overhead = 
    4.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((L->n_sn)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((L->n_sn)*sizeof(int));  /* singlefile matrix arrays */
  printf("real memory overhead %.0f memory %.0f\n",memory_overhead,memory);
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  wtime = taucs_wtime();
  ctime = taucs_ctime();
  
  taucs_io_append(handle,0,1,1,TAUCS_INT,&(L->n_sn));
  taucs_io_append(handle,1,1,L->n_sn+1,TAUCS_INT,L->first_child);
  taucs_io_append(handle,2,1,L->n_sn+1,TAUCS_INT,L->next_child);
  taucs_io_append(handle,3,1,L->n_sn,TAUCS_INT,L->sn_size);
  taucs_io_append(handle,4,1,L->n_sn,TAUCS_INT,L->sn_up_size);
  /*taucs_io_append(handle,5,1,1,TAUCS_INT,&(L->n));*/
  taucs_io_append(handle,6,1,1,TAUCS_INT,&(A->flags));
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Prepare L = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  map  = (int*)taucs_malloc((A->n+1)*sizeof(int));

  sn_in_core = (int*)taucs_malloc((L->n_sn+1)*sizeof(int));
  for(i=0;i<=L->n_sn;i++)
    sn_in_core[i] = 0;
  wtime = taucs_wtime();
  ctime = taucs_ctime();
  if(recursive_compute_supernodes_ll_in_core(L->n_sn,
					     TRUE,
					     (memory - memory_overhead)/3.0,
					     sn_in_core,
					     L)<0.0) {
    ooc_supernodal_factor_free(L);
    taucs_free(map);
    return NULL;
  }
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Scheduling = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  wtime = taucs_wtime();
  ctime = taucs_ctime();

  if (recursive_leftlooking_supernodal_factor_llt_ooc((L->n_sn),  
						      TRUE, 
						      map,
						      sn_in_core,
						      handle,
						      A,L)) {
    ooc_supernodal_factor_free(L);
    taucs_free(map);
    return NULL;
  }

  taucs_printf("\t\tSupernodal Left-Looking INFO: nreads = %.3f bytes_read = %.3f time_read = %10.3f seconds \n",
	       handle->nreads,handle->bytes_read,handle->read_time); 

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSupernodal Left-Looking LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_free(map);
  taucs_free(sn_in_core);
  taucs_io_close(handle);
  ooc_supernodal_factor_free(L);

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Cleanup = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  return (void*) "/tmp/taucs.L";
}
#endif

/******************** OOC SOLVE **************************/

static void 
recursive_supernodal_solve_l_ooc(int sn,       /* this supernode */
				 int is_root,  /* is v the root? */
				 taucs_io_handle* handle,
				 int n_sn,
				 int* first_child, int* next_child,
				 int** sn_struct, int* sn_sizes, int* sn_up_sizes,
				 taucs_datatype x[], taucs_datatype b[],
				 taucs_datatype t[])
{
  int child;
  int  sn_size; /* number of rows/columns in the supernode    */
  int  up_size; /* number of rows that this supernode updates */
  int    ione = 1;
  /*
  double done = 1.0;
  double dzero = 0.0;
  */
  taucs_datatype* xdense;
  taucs_datatype* bdense;
  double  flops;
  int i,j,ip,jp;
  taucs_datatype* sn_block;
  taucs_datatype* up_block = NULL; /* warning */

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    recursive_supernodal_solve_l_ooc(child,
				     FALSE,
				     handle,
				     n_sn,
				     first_child,next_child,
				     sn_struct,sn_sizes,sn_up_sizes,
				     x,b,t);
  }

  if(!is_root) {
    sn_size = sn_sizes[sn];
    up_size = sn_up_sizes[sn] - sn_sizes[sn];

    sn_struct[sn] = (int*)taucs_malloc((sn_size+up_size)*sizeof(int));
    taucs_io_read(handle,IO_BASE+sn,1,sn_size+up_size,TAUCS_INT,sn_struct[sn]);
    
    sn_block = (taucs_datatype*)taucs_calloc(sn_size*sn_size,sizeof(taucs_datatype));
    taucs_io_read(handle,IO_BASE+n_sn+2*sn,
		  sn_size,
		  sn_size ,
		  TAUCS_CORE_DATATYPE,sn_block);
   
    if (up_size > 0 && sn_size > 0) {
      up_block = (taucs_datatype*)taucs_calloc(up_size*sn_size,sizeof(taucs_datatype));
      taucs_io_read(handle,IO_BASE+n_sn+2*sn+1,
		    up_size,
		    sn_size ,
		    TAUCS_CORE_DATATYPE,up_block);
    }

    flops = ((double)sn_size)*((double)sn_size) 
      + 2.0*((double)sn_size)*((double)up_size);

    if (flops > BLAS_FLOPS_CUTOFF) {
      xdense = t;
      bdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	xdense[i] = b[ sn_struct[ sn ][ i ] ];
      for (i=0; i<up_size; i++)
	bdense[i] = taucs_zero;
      
      taucs_trsm ("Left",
	     "Lower",
	     "No Conjugate",
	     "No unit diagonal",
	     &sn_size,&ione,
	     &taucs_one_const,
	     sn_block,&sn_size,
	     xdense , &sn_size);
      
      if (up_size > 0 && sn_size > 0) 
	taucs_gemm ("No Conjugate","No Conjugate",
		    &up_size, &ione, &sn_size,
		    &taucs_one_const,
		    up_block,&up_size,
		    xdense       ,&sn_size,
		    &taucs_zero_const,
		    bdense       ,&up_size);
      
      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = xdense[i];
      for (i=0; i<up_size; i++)
	/*b[ sn_struct[ sn ][ sn_size + i ] ] -= bdense[i];*/
	b[ sn_struct[ sn ][ sn_size + i ] ] =
	  taucs_sub(b[ sn_struct[ sn ][ sn_size + i ] ] , bdense[i]);

    } else if (sn_size > SOLVE_DENSE_CUTOFF) {

      xdense = t;
      bdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	xdense[i] = b[ sn_struct[ sn ][ i ] ];
      for (i=0; i<up_size; i++)
	bdense[i] = taucs_zero;
      
      for (jp=0; jp<sn_size; jp++) {
	/*xdense[jp] = xdense[jp] / sn_block[ sn_size*jp + jp];*/
	xdense[jp] = taucs_div(xdense[jp] , sn_block[ sn_size*jp + jp]);

	for (ip=jp+1; ip<sn_size; ip++) {
	  /*xdense[ip] -= xdense[jp] * sn_block[ sn_size*jp + ip];*/
	  xdense[ip] = taucs_sub(xdense[i],
				 taucs_mul(xdense[jp] , sn_block[ sn_size*jp + ip]));
	}
      }

      for (jp=0; jp<sn_size; jp++) {
	for (ip=0; ip<up_size; ip++) {
	  /*bdense[ip] += xdense[jp] * up_block[ up_size*jp + ip];*/
	  bdense[ip] = taucs_add(bdense[ip],
				 taucs_mul(xdense[jp] , up_block[ up_size*jp + ip]));
	}
      }

      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = xdense[i];
      for (i=0; i<up_size; i++)
	/*b[ sn_struct[ sn ][ sn_size + i ] ] -= bdense[i];*/
	b[ sn_struct[ sn ][ sn_size + i ] ] =
	  taucs_sub(b[ sn_struct[ sn ][ sn_size + i ] ] , bdense[i]);
      
    } else {
      for (jp=0; jp<sn_size; jp++) {
	j = sn_struct[sn][jp];
	/*x[j] = b[j] / sn_block[ sn_size*jp + jp];*/
	x[j] = taucs_div(b[j] , sn_block[ sn_size*jp + jp]);
	for (ip=jp+1; ip<sn_size; ip++) {
	  i = sn_struct[sn][ip];
	  /*b[i] -= x[j] * sn_block[ sn_size*jp + ip];*/
	  b[i] = taucs_sub(b[i],
			   taucs_mul(x[j] , sn_block[ sn_size*jp + ip]));
	}

	for (ip=0; ip<up_size; ip++) {
	  i = sn_struct[sn][sn_size + ip];
	  /*b[i] -= x[j] * up_block[ up_size*jp + ip];*/
	  b[i] = taucs_sub(b[i],
			   taucs_mul(x[j] , up_block[ up_size*jp + ip]));
	}
      }
    }
    taucs_free(sn_struct[sn]);
    taucs_free(sn_block);
    if (up_size > 0 && sn_size > 0) taucs_free(up_block);
    sn_struct[sn] = NULL;
    sn_block = NULL;
    up_block = NULL;
  }
}

static void 
recursive_supernodal_solve_lt_ooc(int sn,       /* this supernode */
				  int is_root,  /* is v the root? */
				  taucs_io_handle* handle,
				  int n_sn,
				  int* first_child, int* next_child,
				  int** sn_struct, int* sn_sizes, int* sn_up_sizes,
				  taucs_datatype x[], taucs_datatype b[],
				  taucs_datatype t[])
{
  int child;
  int  sn_size; /* number of rows/columns in the supernode    */
  int  up_size; /* number of rows that this supernode updates */
  int    ione = 1;
  taucs_datatype* xdense;
  taucs_datatype* bdense;
  double  flops;
  int i,j,ip,jp;
  taucs_datatype* sn_block;
  taucs_datatype* up_block = NULL; /* warning */

  if(!is_root) {

    sn_size = sn_sizes[sn];
    up_size = sn_up_sizes[sn]-sn_sizes[sn];

    sn_struct[sn] = (int*)taucs_malloc((sn_size+up_size)*sizeof(int));
    taucs_io_read(handle,IO_BASE+sn,1,sn_size+up_size,TAUCS_INT,sn_struct[sn]);
    
    sn_block = (taucs_datatype*)taucs_calloc(sn_size*sn_size,sizeof(taucs_datatype));
    taucs_io_read(handle,IO_BASE+n_sn+2*sn,
		  sn_size,
		  sn_size ,
		  TAUCS_CORE_DATATYPE,sn_block);
    if (up_size > 0 && sn_size > 0){
      up_block = (taucs_datatype*)taucs_calloc(up_size*sn_size,sizeof(taucs_datatype));
      taucs_io_read(handle,IO_BASE+n_sn+2*sn+1,
		    up_size,
		    sn_size ,
		    TAUCS_CORE_DATATYPE,up_block);
    }

    flops = ((double)sn_size)*((double)sn_size) 
      + 2.0*((double)sn_size)*((double)up_size);

    if (flops > BLAS_FLOPS_CUTOFF) {
      bdense = t;
      xdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	bdense[i] = b[ sn_struct[ sn][ i ] ];
      for (i=0; i<up_size; i++)
	xdense[i] = x[ sn_struct[sn][sn_size+i] ];
      
      if (up_size > 0 && sn_size > 0)
	taucs_gemm ("Conjugate","No Conjugate",
	       &sn_size, &ione, &up_size,
	       &taucs_minusone_const,
	       up_block,&up_size,
	       xdense       ,&up_size,
	       &taucs_one_const,
	       bdense       ,&sn_size);
            
      taucs_trsm ("Left",
	     "Lower",
	     "Conjugate",
	     "No unit diagonal",
	     &sn_size,&ione,
	     &taucs_one_const,
	     sn_block,&sn_size,
	     bdense       ,&sn_size);
      
      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = bdense[i];
    
    } else if (sn_size > SOLVE_DENSE_CUTOFF) {
      bdense = t;
      xdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	bdense[i] = b[ sn_struct[ sn][ i ] ];
      
      for (i=0; i<up_size; i++)
	xdense[i] = x[ sn_struct[sn][sn_size+i] ];
     
      for (ip=sn_size-1; ip>=0; ip--) {
	for (jp=0; jp<up_size; jp++) {
	  /*bdense[ip] -= xdense[jp] * up_block[ up_size*ip + jp];*/
	  bdense[ip] = taucs_sub(bdense[ip],
				 taucs_mul(xdense[jp] , up_block[ up_size*ip + jp]));
	}
      }
     
      for (ip=sn_size-1; ip>=0; ip--) {
	for (jp=sn_size-1; jp>ip; jp--) {
	  /*bdense[ip] -= bdense[jp] * sn_block[ sn_size*ip + jp];*/
	  bdense[ip] = taucs_sub(bdense[ip],
				 taucs_mul(bdense[jp] , sn_block[ sn_size*ip + jp]));
	}
	/*bdense[ip] = bdense[ip] / sn_block[ sn_size*ip + ip];*/
	bdense[ip] = taucs_div(bdense[ip] , sn_block[ sn_size*ip + ip]);
      }

      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = bdense[i];
    
    } else {
      for (ip=sn_size-1; ip>=0; ip--) {
	i = sn_struct[sn][ip];

	for (jp=0; jp<up_size; jp++) {
	  j = sn_struct[sn][sn_size + jp];
	  /*b[i] -= x[j] * up_block[ up_size*ip + jp];*/
	  b[i] = taucs_sub(b[i],
			   taucs_mul(x[j] , up_block[ up_size*ip + jp]));
	}

	for (jp=sn_size-1; jp>ip; jp--) {
	  j = sn_struct[sn][jp];
	  /*b[i] -= x[j] * sn_block[ sn_size*ip + jp];*/
	  b[i] = taucs_sub(b[i],
			   taucs_mul(x[j] , sn_block[ sn_size*ip + jp]));
	}
	/*x[i] = b[i] / sn_block[ sn_size*ip + ip];*/
	x[i] = taucs_div(b[i] , sn_block[ sn_size*ip + ip]);
      }

    }
    taucs_free(sn_struct[sn]);
    taucs_free(sn_block);
    if (up_size > 0 && sn_size > 0) taucs_free(up_block);
    sn_struct[sn] = NULL;
    sn_block = NULL;
    up_block = NULL;
  }

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    recursive_supernodal_solve_lt_ooc(child,
				      FALSE,
				      handle,
				      n_sn,
				      first_child,next_child,
				      sn_struct,sn_sizes,sn_up_sizes,
				      x,b,t);
  }
}

int taucs_dtl(ooc_solve_llt)(void* vL,
			     void* vx, void* vb)
{
  taucs_io_handle* handle = (taucs_io_handle*) vL;
  taucs_datatype* x = (taucs_datatype*) vx;
  taucs_datatype* b = (taucs_datatype*) vb;
  /*  char* filename = (char*) vL; */
  supernodal_factor_matrix* L;

  taucs_datatype* y;
  taucs_datatype* t; /* temporary vector */
  int     i;

  L = multifrontal_supernodal_create();
  /* READ n, n_sn, first_child, next_child*/
  /*
  if(!MULTIFILE)
    handle = taucs_io_open_singlefile(filename);
  else
    handle = taucs_io_open_multifile(filename);
  */
  taucs_io_read(handle,5,1,1,TAUCS_INT,&(L->n));
  taucs_io_read(handle,0,1,1,TAUCS_INT,&(L->n_sn));
  L->sn_struct = (int**)taucs_malloc((L->n_sn  )*sizeof(int*));
  L->sn_blocks = (taucs_datatype**)taucs_malloc((L->n_sn  )*sizeof(taucs_datatype*));
  L->up_blocks = (taucs_datatype**)taucs_malloc((L->n_sn  )*sizeof(taucs_datatype*));
  L->sn_size   = (int*) taucs_malloc((L->n_sn+1)*sizeof(int));
  L->sn_up_size   = (int*) taucs_malloc((L->n_sn+1)*sizeof(int));
  L->first_child = (int*) taucs_malloc((L->n_sn+1)*sizeof(int));
  L->next_child  = (int*) taucs_malloc((L->n_sn+1)*sizeof(int));
  taucs_io_read(handle,1,1,L->n_sn+1,TAUCS_INT,L->first_child);
  taucs_io_read(handle,2,1,L->n_sn+1,TAUCS_INT,L->next_child);
  taucs_io_read(handle,3,1,L->n_sn,TAUCS_INT,L->sn_size);
  taucs_io_read(handle,4,1,L->n_sn,TAUCS_INT,L->sn_up_size);
  /*  for (i=0; i<L->n_sn; i++) {
    L->sn_struct[i] = (int*)taucs_malloc(L->sn_up_size[i]*sizeof(int));
      taucs_io_read(handle,IO_BASE+i,1,L->sn_up_size[i],TAUCS_INT,L->sn_struct[i]);
      }*/
  
  for(i=0;i<L->n_sn;i++){
    L->sn_struct[i] = NULL;
    L->sn_blocks[i] = NULL;
    L->up_blocks[i] = NULL;
  }
  
  y = (taucs_datatype*) taucs_malloc((L->n) * sizeof(taucs_datatype));
  t = (taucs_datatype*) taucs_malloc((L->n) * sizeof(taucs_datatype));
  if (!y || !t) {
    taucs_free(y);
    taucs_free(t);
    taucs_printf("leftlooking_supernodal_solve_llt: out of memory\n");
    return -1;
  }

  for (i=0; i<L->n; i++) x[i] = b[i];

  recursive_supernodal_solve_l_ooc (L->n_sn,
				    TRUE,  /* this is the root */
				    handle,
				    L->n_sn,
				    L->first_child, L->next_child,
				    L->sn_struct,L->sn_size,L->sn_up_size,
				    y, x, t);

  recursive_supernodal_solve_lt_ooc(L->n_sn,
				    TRUE,  /* this is the root */
				    handle,
				    L->n_sn,
				    L->first_child, L->next_child,
				    L->sn_struct,L->sn_size,L->sn_up_size,
				    x, y, t);

  taucs_free(y);
  taucs_free(t);
  ooc_supernodal_factor_free(L);
  return 0;
}
/*******************************************************************/
/**                     OOC Panelize Factor                       **/
/*******************************************************************/

static double
recursive_smart_panelize_ooc_supernodes(int sn,       /* this supernode */
					int is_root,  /* is v the root? */
					double global_mem,
					int* curr_panel,
					int* sn_in_core,
					int* sn_to_panel_map,
					supernodal_factor_matrix* L)
{
  int  child;
  double this_sn_mem = 0.0;
  double max_child_mem = 0.0;
  double curr_child_mem = 0.0;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;


  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if(!sn_in_core[child]){
      curr_child_mem = recursive_smart_panelize_ooc_supernodes(child,
							       FALSE,
							       global_mem,
							       curr_panel,
							       sn_in_core,
							       sn_to_panel_map,
							       L);
      if(curr_child_mem > max_child_mem) max_child_mem = curr_child_mem;
    }
  }

  if (!is_root){
    this_sn_mem = 
      (double) (L->sn_size)[sn] * (double) (L->sn_up_size)[sn] * (double) sizeof(taucs_datatype) 
      + (double) (L->sn_up_size)[sn] * (double) sizeof(int);
    if(max_child_mem+this_sn_mem < global_mem){
      sn_to_panel_map[sn] = *curr_panel;
      return (max_child_mem+this_sn_mem);
    } else {
      (*curr_panel)++;
      sn_to_panel_map[sn] = *curr_panel;
      return this_sn_mem;
    }
  } 

  /* reached only at the root */
  return 0.0;
}

static double
recursive_dumb_panelize_ooc_supernodes(int sn,       /* this supernode */
				       int is_root,  /* is v the root? */
				       int* curr_panel,
				       int* sn_in_core,
				       int* sn_to_panel_map,
				       supernodal_factor_matrix* L)
{
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if(!sn_in_core[child]) 
      (void) recursive_dumb_panelize_ooc_supernodes(child,
						    FALSE,
						    curr_panel,
						    sn_in_core,
						    sn_to_panel_map,
						    L);
  }

  if (!is_root){
    (*curr_panel)++;
    sn_to_panel_map[sn]=*curr_panel;
  }

  return 0.0;
}

static double
recursive_panelize_ooc_supernodes(int sn,       /* this supernode */
				  int is_root,  /* is v the root? */
				  double global_mem,
				  double avail_mem,
				  int* curr_panel,
				  int* sn_in_core,
				  int* sn_to_panel_map,
				  supernodal_factor_matrix* L)
{
  int  child;
  double curr_avail_mem = avail_mem;
  double this_sn_mem = 0.0;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;


  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if(!sn_in_core[child]) 
      curr_avail_mem=recursive_panelize_ooc_supernodes(child,
						       FALSE,
						       global_mem,
						       curr_avail_mem,
						       curr_panel,
						       sn_in_core,
						       sn_to_panel_map,
						       L);
  }

  if (!is_root){
    this_sn_mem = (double)(L->sn_size)[sn]*(double)(L->sn_up_size)[sn]*sizeof(taucs_datatype)+(double)(L->sn_up_size)[sn]*sizeof(int);
    if(curr_avail_mem-this_sn_mem>0.0){
      sn_to_panel_map[sn]=*curr_panel;
      curr_avail_mem -= this_sn_mem;
    } else {
      (*curr_panel)++;
      sn_to_panel_map[sn]=*curr_panel;
      curr_avail_mem = global_mem-this_sn_mem;
    }
  }
  return curr_avail_mem;
}

static void
recursive_leftlooking_supernodal_update_panel_ooc(int J,int K,
						  int bitmap[],
						  int* sn_to_panel_map,
						  taucs_datatype* dense_update_matrix,
						  taucs_io_handle* handle,
						  taucs_ccs_matrix* A,
						  supernodal_factor_matrix* L)
{
  int i,j,ir,ii;
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  int sn_size_child = (L->sn_size)[K];
  int sn_up_size_child = (L->sn_up_size)[K];
  int exist_upd=0;
  int first_row=0;
  int row_count=0;
  int updated_panel = sn_to_panel_map[J];
  int curr_updated_sn = J;
  int PK,M,N,LDA,LDB,LDC;
 
  if(L->sn_up_size[K]-L->sn_size[K]>0){
    
    if(!(L->sn_struct)[K]){
      L->sn_struct[K] = (int*)taucs_malloc(sn_up_size_child*sizeof(int));
      taucs_io_read(handle,IO_BASE+K,1,sn_up_size_child,TAUCS_INT,L->sn_struct[K]);
    }

    for(i=sn_size_child;i<sn_up_size_child;i++){
      /* We want to update only supernodes great than J with using of
	 sort of indices */
      if(L->col_to_sn_map[L->sn_struct[K][i]]<curr_updated_sn)
	continue;

      /* is this row index included in the columns of curr_updated_sn ? */
      if(L->col_to_sn_map[L->sn_struct[K][i]]==curr_updated_sn){
	if(!exist_upd) first_row = i;
	row_count++;
	exist_upd = 1;
      }
     
      /* if end of update to curr_updated_sn or edge condition */
      if(L->col_to_sn_map[L->sn_struct[K][i]]!=curr_updated_sn||
	 i==sn_up_size_child-1)
	if(exist_upd){ 
	  if(!(L->up_blocks)[K]){
	    (L->up_blocks)[K] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[K]-(L->sn_size)[K])
							*((L->sn_size)[K]),
							sizeof(taucs_datatype));
	    taucs_io_read(handle,IO_BASE+L->n_sn+2*K+1,
			  (L->sn_up_size)[K]-(L->sn_size)[K],
			  (L->sn_size)[K] ,
			  TAUCS_CORE_DATATYPE,(L->up_blocks)[K]);
	  }
	  LDA = LDB = (L->sn_up_size)[K]-(L->sn_size)[K];
	  M  = sn_up_size_child - first_row ;
	  LDC =  M;
	  N  = row_count; 
	  PK = L->sn_size[K];    
	  
	  /* The GEMM code computes on the upper triangle of the trapezoidal
	     matrix, which is junk. */
	  /*
             taucs_gemm ("No Conjugate",
		         "Conjugate",
		         &M,&N,&PK,
		         &taucs_one_const,
		         &(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		         &(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		         &taucs_zero_const,
		         dense_update_matrix,&LDC);
	  */

	  /* This is the HERK+GEMM fix by Elad */
	  taucs_herk ("Lower",
		      "No Conjugate",
		      &N,&PK,
		      &taucs_one_real_const,
		      &(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		      &taucs_zero_real_const,
		      dense_update_matrix,&LDC);

	  if(M-N > 0)
	    {
	      int newM = M - N;
	      
	      taucs_gemm ("No Conjugate",
			  "Conjugate",
			  &newM,&N,&PK,
			  &taucs_one_const,
			  &(L->up_blocks[K][first_row-sn_size_child+N]),&LDA,
			  &(L->up_blocks[K][first_row-sn_size_child]),&LDB,
			  &taucs_zero_const,
			  dense_update_matrix+N,&LDC);
	    }
	  /* end of GEMM/HERK+GEMM fix */ 

	  if(!(L->sn_blocks)[curr_updated_sn])
	    (L->sn_blocks)[curr_updated_sn] = 
	      (taucs_datatype*)taucs_calloc(((L->sn_size)[curr_updated_sn])
				      *((L->sn_size)[curr_updated_sn]),
				      sizeof(taucs_datatype));
    
	  if(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn]>0)
	    if(!(L->up_blocks)[curr_updated_sn]) 
	      (L->up_blocks)[curr_updated_sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[curr_updated_sn]-(L->sn_size)[curr_updated_sn]) *((L->sn_size)[curr_updated_sn]),sizeof(taucs_datatype));

	  if(!L->sn_struct[curr_updated_sn]){
	    L->sn_struct[curr_updated_sn] = (int*)taucs_malloc((L->sn_up_size)[curr_updated_sn]*sizeof(int));
	    taucs_io_read(handle,
			  IO_BASE+curr_updated_sn,
			  1,(L->sn_up_size)[curr_updated_sn],
			  TAUCS_INT,L->sn_struct[curr_updated_sn]);
	  }

	  for(ii=0;ii<L->sn_size[curr_updated_sn];ii++) {
	    bitmap[L->sn_struct[curr_updated_sn][ii]]=ii+1;
	  }

	  for(ii=L->sn_size[curr_updated_sn];ii<L->sn_up_size[curr_updated_sn];ii++){
	    bitmap[L->sn_struct[curr_updated_sn][ii]] = ii-L->sn_size[curr_updated_sn]+1;
	  }
	   
	  assert((double)row_count*(double)LDC < 2048.0*1024.0*1024.0);
	  for(j=0;j<row_count;j++)
	    for(ir=j;ir<row_count;ir++){
	      /*
		L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] -= dense_update_matrix[j*LDC+ir];*/
	      L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] =
		taucs_sub(L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] , dense_update_matrix[j*LDC+ir]);

	      /* to find overflows */
	      assert((double)(bitmap[L->sn_struct[K][first_row+j]]-1)*(double)L->sn_size[curr_updated_sn] < 2048.0*1024.0*1024.0);
	     }
	  for(j=0;j<row_count;j++)
	    for(ir=row_count;ir<M;ir++){
	      /*L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] -= dense_update_matrix[j*LDC+ir];*/

	      L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] =
		taucs_sub(L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] , dense_update_matrix[j*LDC+ir]);

	      /* to find overflow */
	      assert((double)(bitmap[L->sn_struct[K][first_row+j]]-1)*(double)(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn]) < 2048.0*1024.0*1024.0);
	     }	
	  for(ii=0;ii<L->sn_up_size[curr_updated_sn];ii++)
	    bitmap[L->sn_struct[curr_updated_sn][ii]]=0;

	  exist_upd = 0;
	  row_count = 0;
	}

      /* is this row index included in the columns of sn in the same panel? */
      if(L->col_to_sn_map[L->sn_struct[K][i]]!=curr_updated_sn)
	if(sn_to_panel_map[L->col_to_sn_map[L->sn_struct[K][i]]]==updated_panel){
	  curr_updated_sn = L->col_to_sn_map[L->sn_struct[K][i]];
	  if(!exist_upd) first_row = i;
	  row_count++;
	  exist_upd = 1;
	  if( i==sn_up_size_child-1)
	    if(exist_upd){ 
	      if(!(L->up_blocks)[K]){
		(L->up_blocks)[K] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[K]-(L->sn_size)[K])*((L->sn_size)[K]),sizeof(taucs_datatype));
		taucs_io_read(handle,IO_BASE+L->n_sn+2*K+1,
			      (L->sn_up_size)[K]-(L->sn_size)[K],
			      (L->sn_size)[K] ,
			      TAUCS_CORE_DATATYPE,(L->up_blocks)[K]);
	      }
	      LDA = LDB = (L->sn_up_size)[K]-(L->sn_size)[K];
	      M  = sn_up_size_child - first_row ;
	      LDC =  M;
	      N  = row_count; 
	      PK = L->sn_size[K];    
	      
   
	      /* The GEMM code computes on the upper triangle of the trapezoidal
		 matrix, which is junk. */
	      /*
	      taucs_gemm ("No Conjugate",
		     "Conjugate",
		     &M,&N,&PK,
		     &taucs_one_const,
		     &(L->up_blocks[K][first_row-sn_size_child]),&LDA,
		     &(L->up_blocks[K][first_row-sn_size_child]),&LDB,
		     &taucs_zero_const,
		     dense_update_matrix,&LDC);
	      */

	      /* This is the HERK+GEMM fix by Elad */
	      taucs_herk ("Lower",
			  "No Conjugate",
			  &N,&PK,
			  &taucs_one_real_const,
			  &(L->up_blocks[K][first_row-sn_size_child]),&LDA,
			  &taucs_zero_real_const,
			  dense_update_matrix,&LDC);
	      
	      if(M-N > 0)
		{
		  int newM = M - N;
		  
		  taucs_gemm ("No Conjugate",
			      "Conjugate",
			      &newM,&N,&PK,
			      &taucs_one_const,
			      &(L->up_blocks[K][first_row-sn_size_child+N]),&LDA,
			      &(L->up_blocks[K][first_row-sn_size_child]),&LDB,
			      &taucs_zero_const,
			      dense_update_matrix+N,&LDC);
		}
	      /* end of GEMM/HERK+GEMM fix */ 
	      
	      if(!(L->sn_blocks)[curr_updated_sn])
		(L->sn_blocks)[curr_updated_sn] = 
		  (taucs_datatype*)taucs_calloc(((L->sn_size)[curr_updated_sn])*((L->sn_size)[curr_updated_sn]),sizeof(taucs_datatype));
    
	      if(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn]>0)
		if(!(L->up_blocks)[curr_updated_sn]) 
		  (L->up_blocks)[curr_updated_sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[curr_updated_sn]-(L->sn_size)[curr_updated_sn]) *((L->sn_size)[curr_updated_sn]),sizeof(taucs_datatype));
	      
	      if(!L->sn_struct[curr_updated_sn]){
		L->sn_struct[curr_updated_sn] = (int*)taucs_malloc((L->sn_up_size)[curr_updated_sn]*sizeof(int));
		taucs_io_read(handle,IO_BASE+curr_updated_sn,1,(L->sn_up_size)[curr_updated_sn],TAUCS_INT,L->sn_struct[curr_updated_sn]);
	      }

	      for(ii=0;ii<L->sn_size[curr_updated_sn];ii++) {
		bitmap[L->sn_struct[curr_updated_sn][ii]]=ii+1;
	      }

	      for(ii=L->sn_size[curr_updated_sn];ii<L->sn_up_size[curr_updated_sn];ii++){
	    bitmap[L->sn_struct[curr_updated_sn][ii]] = ii-L->sn_size[curr_updated_sn]+1;
	      }
	  assert((double)row_count*(double)LDC < 2048.0*1024.0*1024.0);
	      for(j=0;j<row_count;j++)
		for(ir=j;ir<row_count;ir++){
		  /*L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] -= dense_update_matrix[j*LDC+ir];*/
		  L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] =
		    taucs_sub(L->sn_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*L->sn_size[curr_updated_sn]+(bitmap[L->sn_struct[K][first_row+ir]]-1)] , dense_update_matrix[j*LDC+ir]);

		  /* for find overflow */
		  assert((double)(bitmap[L->sn_struct[K][first_row+j]]-1)*(double)L->sn_size[curr_updated_sn] < 2048.0*1024.0*1024.0);
		}

	      for(j=0;j<row_count;j++)
		for(ir=row_count;ir<M;ir++){

		  /*L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] -= dense_update_matrix[j*LDC+ir];*/

		  L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] =
		    taucs_sub(L->up_blocks[curr_updated_sn][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn])+(bitmap[L->sn_struct[K][ir+first_row]]-1)] , dense_update_matrix[j*LDC+ir]);

		  /* for find overflow */
		  assert((double)(bitmap[L->sn_struct[K][first_row+j]]-1)*(double)(L->sn_up_size[curr_updated_sn]-L->sn_size[curr_updated_sn]) < 2048.0*1024.0*1024.0);
		}	
	      for(ii=0;ii<L->sn_up_size[curr_updated_sn];ii++)
		bitmap[L->sn_struct[curr_updated_sn][ii]]=0;

	      exist_upd = 0;
	      row_count = 0;

	    }
	}
      
    }
    
  }

  /* free update sn from memory */
  taucs_free((L->up_blocks   )[K]);
  (L->up_blocks   )[K] = NULL;
  taucs_free( L->sn_struct[K]);
  L->sn_struct[K] = NULL;

  if(first_row&&sn_to_panel_map[J]!=sn_to_panel_map[K]){   
  
    for (child = first_child[K]; child != -1; child = next_child[child]) {
      recursive_leftlooking_supernodal_update_panel_ooc(J,child,
							bitmap,
							sn_to_panel_map,
							dense_update_matrix,
							handle,A,L);
    }
  } 

}


static int
recursive_leftlooking_supernodal_factor_panel_llt_ooc
(int sn,    /* this supernode */
 int father_sn,
 int is_root,/* is sn the root?*/
 int* map,
 int* sn_in_core,
 int* sn_to_panel_map,
 int*  panel_max_size,
 taucs_io_handle* handle,
 taucs_ccs_matrix* A,
 supernodal_factor_matrix* L)
{
  int  child;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  taucs_datatype* dense_update_matrix = NULL;
  
  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if(sn_in_core[child]){

      if (recursive_read_L_cols(child,
				FALSE,
				handle,
				L)) {
	return -1;
      }
      if (recursive_leftlooking_supernodal_factor_llt(child,
						      FALSE,
						      map,
						      A,L)) {
	/* failure */
	return -1;
      }
      if (recursive_append_L(child,
			     FALSE,
			     handle,
			     L)) {
	return -1;
      }

    }
    else
      if (recursive_leftlooking_supernodal_factor_panel_llt_ooc(child,
								sn,
								FALSE,
								map,
								sn_in_core,
								sn_to_panel_map,
								panel_max_size,
								handle,
								A,L)) {
	/* failure */
	return -1;
      }
  }

  if (!is_root) { 
    if(!(L->sn_blocks)[sn])
      (L->sn_blocks)[sn] = 
	(taucs_datatype*)taucs_calloc(((L->sn_size)[sn])*((L->sn_size)[sn]),sizeof(taucs_datatype));
    
    if(L->sn_up_size[sn]-L->sn_size[sn]>0)
      if(!(L->up_blocks)[sn]) 
	(L->up_blocks)[sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn]) *((L->sn_size)[sn]),sizeof(taucs_datatype));

    if(!L->sn_struct[sn]){
      L->sn_struct[sn] = (int*)taucs_malloc((L->sn_up_size)[sn]*sizeof(int));
      taucs_io_read(handle,IO_BASE+sn,1,(L->sn_up_size)[sn],TAUCS_INT,L->sn_struct[sn]);
    }

    if (!dense_update_matrix) 
      dense_update_matrix = (taucs_datatype*) taucs_calloc(panel_max_size[sn_to_panel_map[sn]],sizeof(taucs_datatype));

    for (child = first_child[sn]; child != -1; child = next_child[child]) {
      if(sn_to_panel_map[sn]!=sn_to_panel_map[child])
	recursive_leftlooking_supernodal_update_panel_ooc(sn,child,
							  map,
							  sn_to_panel_map,
							  dense_update_matrix,
							  handle,A,L);
    }

    if (leftlooking_supernodal_front_factor(sn,
					    map,
					    A,
					    L)) {
      /* nonpositive pivot */
      return -1;
    }
    
    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn,
		    L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->sn_blocks[sn]);

    taucs_io_append(handle,IO_BASE+L->n_sn+2*sn+1,
		    L->sn_up_size[sn] - L->sn_size[sn],L->sn_size[sn],
		    TAUCS_CORE_DATATYPE,L->up_blocks[sn]);
 
    
    if(sn_to_panel_map[sn]==sn_to_panel_map[father_sn])
      recursive_leftlooking_supernodal_update_panel_ooc(father_sn,sn,
							map,
							sn_to_panel_map,
							dense_update_matrix,
							handle,A,L);
    taucs_free(dense_update_matrix);
    taucs_free((L->sn_blocks)[sn]);
    taucs_free((L->up_blocks)[sn]); 
    taucs_free((L->sn_struct)[sn]);
    (L->sn_blocks   )[sn] = NULL;      
    (L->up_blocks   )[sn] = NULL;
    (L->sn_struct   )[sn] = NULL;
  }

  return 0;
}

int taucs_dtl(ooc_factor_llt)(taucs_ccs_matrix* A, 
			      taucs_io_handle* handle,
			      double memory)
{
  supernodal_factor_matrix* L;
  int i;
  int* map;
  int* sn_in_core;
  int* sn_to_panel_map;
  int* panel_max_size;
  int n_pn=0;
  double wtime, ctime;

  /*
  int j,ip,jp;
  int sn,p;
  int current_index=0;
  */
  double memory_overhead;
  double  max_multiple=0.0;
  int ind_max_mult = 0;

  /* compute fixed memory overhead */
  
  memory_overhead = 
    4.0*(double)((A->n)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((A->n)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((A->n)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((A->n)*sizeof(int));  /* singlefile matrix arrays */

  taucs_printf("\t\tOOC memory overhead bound %.0lf MB (out of %.0lf MB available)\n",
	       memory_overhead/1048576.0,memory/1048576.0);

  taucs_printf(">>> 1\n");

  if ( memory - memory_overhead < 
       2.0*(double)((A->n)*sizeof(taucs_datatype)) + 
       2.0*(double)((A->n)*sizeof(int)) ) {
    taucs_printf("\t\ttaucs_ccs_factor_llt_ll_ooc: not enough memory\n");
    return -1;
    }

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();
  /*
  if(!MULTIFILE)
    handle = taucs_io_create_singlefile("/tmp/taucs.L");
  else
    handle = taucs_io_create_multifile("/tmp/taucs.L");
  */
  taucs_io_append(handle,5,1,1,TAUCS_INT,&(A->n));

  taucs_ccs_ooc_symbolic_elimination(A,L,
				     TRUE /* sort row indices */,
				     TRUE /* return col_to_sn_map */,
				     (memory - memory_overhead)/3.0,
				     ooc_sn_struct_handler,handle);
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  /* we now compute an exact memory overhead bound using n_sn */
  memory_overhead = 
    /*    1.0*(double)((L->n)*sizeof(int)) + */ /* integer vector in L */
    4.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((L->n_sn)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((L->n_sn)*sizeof(int));  /* singlefile matrix arrays */

  taucs_printf("\t\tOOC actual memory overhead %.0lf MB (out of %.0lf MB available)\n",
	       memory_overhead/1048576.0,memory/1048576.0);
 
  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_io_append(handle,0,1,1,TAUCS_INT,&(L->n_sn));
  taucs_io_append(handle,1,1,L->n_sn+1,TAUCS_INT,L->first_child);
  taucs_io_append(handle,2,1,L->n_sn+1,TAUCS_INT,L->next_child);
  taucs_io_append(handle,3,1,L->n_sn,TAUCS_INT,L->sn_size);
  taucs_io_append(handle,4,1,L->n_sn,TAUCS_INT,L->sn_up_size);
  /*taucs_io_append(handle,5,1,1,TAUCS_INT,&(L->n));*/
  taucs_io_append(handle,6,1,1,TAUCS_INT,&(A->flags));

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Prepare L = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  map  = (int*)taucs_malloc((A->n+1)*sizeof(int));
  sn_in_core = (int*)taucs_malloc((L->n_sn+1)*sizeof(int));
  sn_to_panel_map = (int*)taucs_malloc((L->n_sn+1)*sizeof(int));
  for(i=0;i<=L->n_sn;i++){
    sn_in_core[i] = 0;
    sn_to_panel_map[i]=-1;
  }

  for(i=0;i<L->n_sn;i++){
    (L->sn_blocks)[i] = NULL;
    (L->up_blocks)[i] = NULL;
    (L->sn_struct)[i] = NULL;
  }

  wtime = taucs_wtime();
  ctime = taucs_ctime();
  if(recursive_compute_supernodes_ll_in_core(L->n_sn,
					     TRUE,
					     (memory - memory_overhead)/3.0,
					     sn_in_core,
					     L)<0.0) {
    ooc_supernodal_factor_free(L);
    taucs_free(sn_in_core);
    taucs_free(sn_to_panel_map);  
    taucs_free(map);
    return -1;
  }
  
 /*  if(recursive_panelize_ooc_supernodes(L->n_sn,
				       TRUE,
				       (memory - memory_overhead)/3.0,
				       (memory - memory_overhead)/3.0,
				       &n_pn,
				       sn_in_core,
				       sn_to_panel_map,
				       L)<0.0){
    ooc_supernodal_factor_free(L);
    taucs_free(sn_in_core);
    taucs_free(sn_to_panel_map);  
    taucs_free(map);
    return -1;
    }*/
 
  taucs_printf("\t\tOOC Supernodal Left-Looking: panel-is-paged\n",n_pn);

  if(recursive_smart_panelize_ooc_supernodes(L->n_sn,
					     TRUE,
					     (memory - memory_overhead)/3.0,
					     &n_pn,
					     sn_in_core,
					     sn_to_panel_map,
					     L)<0.0){
    ooc_supernodal_factor_free(L);
    taucs_free(sn_in_core);
    taucs_free(sn_to_panel_map);  
    taucs_free(map);
    return -1;
    }

  /* it will be at least one panel even empty */
  n_pn++; 

  /*for(i=0;i<L->n_sn;i++){
    taucs_printf("sn_in_core[%d] = %d\n",i,sn_in_core[i]);
    taucs_printf("sn_to_panel_map[%d] = %d\n",i,sn_to_panel_map[i]);
    }*/
  taucs_printf("\t\tOOC Supernodal Left-Looking: %d panels\n",n_pn);
  /* compute max dense matrix size for every panel */
  panel_max_size = (int*)taucs_calloc(n_pn,sizeof(int));
  for(i=0;i<L->n_sn;i++){
    if((double)L->sn_up_size[i]*(double)L->sn_size[i]>max_multiple){
      max_multiple = (double)L->sn_up_size[i]*(double)L->sn_size[i]; 
      ind_max_mult = i;
      }
    if(sn_to_panel_map[i]!=-1){
      if(L->sn_up_size[i]*L->sn_size[i]>panel_max_size[sn_to_panel_map[i]]) 
	panel_max_size[sn_to_panel_map[i]] = L->sn_up_size[i]*L->sn_size[i];
    }
  }
  /*
  taucs_printf("debug***: L->n_sn = %d max(sn_size*sn_up_size) = %lf sn_size[%d] = %d sn_up_size[%d] = %d\n ",L->n_sn,max_multiple,ind_max_mult,L->sn_size[ind_max_mult],ind_max_mult,L->sn_up_size[ind_max_mult]);
  */

  /*  for(i=0;i<n_pn;i++)
      taucs_printf(" panel_max_size[%d] = %d\n",i, panel_max_size[i]);*/
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Scheduling = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  wtime = taucs_wtime();
  ctime = taucs_ctime();


  if (recursive_leftlooking_supernodal_factor_panel_llt_ooc(L->n_sn,
							    L->n_sn,  
							    TRUE, 
							    map,
							    sn_in_core,
							    sn_to_panel_map,
							    panel_max_size,
							    handle,
							    A,L)) {
    ooc_supernodal_factor_free(L);
    taucs_free(map);
    return -1;
  }
 
 taucs_printf("\t\tOOC Supernodal Left-Looking:\n");
 taucs_printf("\t\t\tread count           = %.0f \n",handle->nreads);
 taucs_printf("\t\t\tread volume (bytes)  = %.2e \n",handle->bytes_read);
 taucs_printf("\t\t\tread time (seconds)  = %.0f \n",handle->read_time);
 taucs_printf("\t\t\twrite count          = %.0f \n",handle->nwrites);
 taucs_printf("\t\t\twrite volume (bytes) = %.2e \n",handle->bytes_written);
 taucs_printf("\t\t\twrite time (seconds) = %.0f \n",handle->write_time);

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_free(map);
  taucs_free(sn_in_core);
  taucs_free(sn_to_panel_map);  
  /*taucs_io_close(handle);*/
  ooc_supernodal_factor_free(L);
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Cleanup = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  /*return (void*) "/tmp/taucs.L";*/
  /*return (void*) "/rabani2/queen/taucs.L";*/

  return 0;
}

/*************************************************************/
/* SAME ROUTINE, WITH CHOICE OF PANELIZATION FOR TESTING     */
/*************************************************************/

int taucs_dtl(ooc_factor_llt_panelchoice)(taucs_ccs_matrix* A, 
					  taucs_io_handle* handle,
					  double memory,
					  int panelization_method)
{
  supernodal_factor_matrix* L;
  int i;
  int* map;
  int* sn_in_core;
  int* sn_to_panel_map;
  int* panel_max_size;
  int n_pn=0;
  double wtime, ctime;
  /*
  int j,ip,jp;
  int sn,p;
  int current_index=0;
  */
  double memory_overhead;
  double  max_multiple=0.0;
  int ind_max_mult = 0;

  /* compute fixed memory overhead */
  
  memory_overhead = 
    4.0*(double)((A->n)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((A->n)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((A->n)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((A->n)*sizeof(int));  /* singlefile matrix arrays */

  taucs_printf("\t\tOOC memory overhead bound %.0lf MB (out of %.0lf MB available)\n",
	       memory_overhead/1048576.0,memory/1048576.0);

  taucs_printf("*** 1\n");

  if ( memory - memory_overhead < 
       2.0*(double)((A->n)*sizeof(taucs_datatype)) + 
       2.0*(double)((A->n)*sizeof(int)) ) {
    taucs_printf("\t\ttaucs_ccs_factor_llt_ll_ooc: not enough memory\n");
    return -1;
    }

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_printf("*** 2\n");

  L = multifrontal_supernodal_create();
  /*
  if(!MULTIFILE)
    handle = taucs_io_create_singlefile("/tmp/taucs.L");
  else
    handle = taucs_io_create_multifile("/tmp/taucs.L");
  */
  taucs_io_append(handle,5,1,1,TAUCS_INT,&(A->n));

  taucs_printf("*** 3\n");

  taucs_ccs_ooc_symbolic_elimination(A,L,
				     TRUE /* sort row indices */,
				     TRUE /* return col_to_sn_map */,
				     (memory - memory_overhead)/3.0,
				     ooc_sn_struct_handler,handle);
  
  taucs_printf("*** 4\n");

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  /* we now compute an exact memory overhead bound using n_sn */
  memory_overhead = 
    /*    1.0*(double)((L->n)*sizeof(int)) + */ /* integer vector in L */
    4.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in L */
    3.0*(double)((L->n_sn)*sizeof(int)) + /* pointer arrays in L  */
    2.0*(double)((L->n_sn)*sizeof(int)) + /* integer vectors in program  */
    4.0*3.0*(double)((L->n_sn)*sizeof(int));  /* singlefile matrix arrays */

  taucs_printf("\t\tOOC actual memory overhead %.0lf MB (out of %.0lf MB available)\n",
	       memory_overhead/1048576.0,memory/1048576.0);
 
  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_io_append(handle,0,1,1,TAUCS_INT,&(L->n_sn));
  taucs_io_append(handle,1,1,L->n_sn+1,TAUCS_INT,L->first_child);
  taucs_io_append(handle,2,1,L->n_sn+1,TAUCS_INT,L->next_child);
  taucs_io_append(handle,3,1,L->n_sn,TAUCS_INT,L->sn_size);
  taucs_io_append(handle,4,1,L->n_sn,TAUCS_INT,L->sn_up_size);
  /*taucs_io_append(handle,5,1,1,TAUCS_INT,&(L->n));*/
  taucs_io_append(handle,6,1,1,TAUCS_INT,&(A->flags));

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Prepare L = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  map  = (int*)taucs_malloc((A->n+1)*sizeof(int));
  sn_in_core = (int*)taucs_malloc((L->n_sn+1)*sizeof(int));
  sn_to_panel_map = (int*)taucs_malloc((L->n_sn+1)*sizeof(int));
  for(i=0;i<=L->n_sn;i++){
    sn_in_core[i] = 0;
    sn_to_panel_map[i]=-1;
  }

  for(i=0;i<L->n_sn;i++){
    (L->sn_blocks)[i] = NULL;
    (L->up_blocks)[i] = NULL;
    (L->sn_struct)[i] = NULL;
  }

  wtime = taucs_wtime();
  ctime = taucs_ctime();
  if(recursive_compute_supernodes_ll_in_core(L->n_sn,
					     TRUE,
					     (memory - memory_overhead)/3.0,
					     sn_in_core,
					     L)<0.0) {
    ooc_supernodal_factor_free(L);
    taucs_free(sn_in_core);
    taucs_free(sn_to_panel_map);  
    taucs_free(map);
    return -1;
  }

  if (panelization_method == 1) {
    taucs_printf("\t\tOOC Supernodal Left-Looking: panel-in-memory\n",n_pn);
    if(recursive_panelize_ooc_supernodes(L->n_sn,
					 TRUE,
					 (memory - memory_overhead)/3.0,
					 (memory - memory_overhead)/3.0,
					 &n_pn,
					 sn_in_core,
					 sn_to_panel_map,
					 L)<0.0){
      ooc_supernodal_factor_free(L);
      taucs_free(sn_in_core);
      taucs_free(sn_to_panel_map);  
      taucs_free(map);
      return -1;
    }
  } 
 
  if (panelization_method == 0) {
    taucs_printf("\t\tOOC Supernodal Left-Looking: panel-is-paged\n",n_pn);
    if(recursive_smart_panelize_ooc_supernodes(L->n_sn,
					       TRUE,
					       (memory - memory_overhead)/3.0,
					       &n_pn,
					       sn_in_core,
					       sn_to_panel_map,
					       L)<0.0){
      ooc_supernodal_factor_free(L);
      taucs_free(sn_in_core);
      taucs_free(sn_to_panel_map);  
      taucs_free(map);
      return -1;
    }
  }

  if (panelization_method == 2) {
    taucs_printf("\t\tOOC Supernodal Left-Looking: panel-is-supernode\n",n_pn);
    if (recursive_dumb_panelize_ooc_supernodes(L->n_sn,
					       TRUE,
					       &n_pn,
					       sn_in_core,
					       sn_to_panel_map,
					       L)<0.0){
      ooc_supernodal_factor_free(L);
      taucs_free(sn_in_core);
      taucs_free(sn_to_panel_map);  
      taucs_free(map);
      return -1;
    }
  }

  /* it will be at least one panel even empty */
  n_pn++; 

  /*for(i=0;i<L->n_sn;i++){
    taucs_printf("sn_in_core[%d] = %d\n",i,sn_in_core[i]);
    taucs_printf("sn_to_panel_map[%d] = %d\n",i,sn_to_panel_map[i]);
    }*/
  taucs_printf("\t\tOOC Supernodal Left-Looking: %d panels\n",n_pn);
  /* compute max dense matrix size for every panel */
  panel_max_size = (int*)taucs_calloc(n_pn,sizeof(int));
  for(i=0;i<L->n_sn;i++){
    if((double)L->sn_up_size[i]*(double)L->sn_size[i]>max_multiple){
      max_multiple = (double)L->sn_up_size[i]*(double)L->sn_size[i]; 
      ind_max_mult = i;
      }
    if(sn_to_panel_map[i]!=-1){
      if(L->sn_up_size[i]*L->sn_size[i]>panel_max_size[sn_to_panel_map[i]]) 
	panel_max_size[sn_to_panel_map[i]] = L->sn_up_size[i]*L->sn_size[i];
    }
  }
  /*
  taucs_printf("debug***: L->n_sn = %d max(sn_size*sn_up_size) = %lf sn_size[%d] = %d sn_up_size[%d] = %d\n ",L->n_sn,max_multiple,ind_max_mult,L->sn_size[ind_max_mult],ind_max_mult,L->sn_up_size[ind_max_mult]);
  */

  /*  for(i=0;i<n_pn;i++)
      taucs_printf(" panel_max_size[%d] = %d\n",i, panel_max_size[i]);*/
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Scheduling = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  wtime = taucs_wtime();
  ctime = taucs_ctime();


  if (recursive_leftlooking_supernodal_factor_panel_llt_ooc(L->n_sn,
							    L->n_sn,  
							    TRUE, 
							    map,
							    sn_in_core,
							    sn_to_panel_map,
							    panel_max_size,
							    handle,
							    A,L)) {
    ooc_supernodal_factor_free(L);
    taucs_free(map);
    return -1;
  }
 
 taucs_printf("\t\tOOC Supernodal Left-Looking:\n");
 taucs_printf("\t\t\tread count           = %.0f \n",handle->nreads);
 taucs_printf("\t\t\tread volume (bytes)  = %.2e \n",handle->bytes_read);
 taucs_printf("\t\t\tread time (seconds)  = %.0f \n",handle->read_time);
 taucs_printf("\t\t\twrite count          = %.0f \n",handle->nwrites);
 taucs_printf("\t\t\twrite volume (bytes) = %.2e \n",handle->bytes_written);
 taucs_printf("\t\t\twrite time (seconds) = %.0f \n",handle->write_time);

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  taucs_free(map);
  taucs_free(sn_in_core);
  taucs_free(sn_to_panel_map);  
  /*taucs_io_close(handle);*/
  ooc_supernodal_factor_free(L);
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tOOC Supernodal Left-Looking Cleanup = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  
  /*return (void*) "/tmp/taucs.L";*/
  /*return (void*) "/rabani2/queen/taucs.L";*/

  return 0;
}

#else /* TAUCS_CORE_GENRAL */

/*************************************************************/
/* generic interfaces to user-callable routines              */
/*************************************************************/

int taucs_ooc_factor_llt(taucs_ccs_matrix* A,
			 taucs_io_handle*  L,
			 double memory)
{
#ifdef TAUCS_CONFIG_DREAL
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dooc_factor_llt(A,L,memory);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (A->flags & TAUCS_SINGLE)
    return taucs_sooc_factor_llt(A,L,memory);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zooc_factor_llt(A,L,memory);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cooc_factor_llt(A,L,memory);
#endif

  assert(0);
  return -1;
}

int taucs_ooc_factor_llt_panelchoice(taucs_ccs_matrix* A,
				     taucs_io_handle*  L,
				     double memory,
				     int panelchoice)
{
#ifdef TAUCS_CONFIG_DREAL
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dooc_factor_llt_panelchoice(A,L,memory,panelchoice);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (A->flags & TAUCS_SINGLE)
    return taucs_sooc_factor_llt_panelchoice(A,L,memory,panelchoice);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zooc_factor_llt_panelchoice(A,L,memory,panelchoice);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cooc_factor_llt_panelchoice(A,L,memory,panelchoice);
#endif

  assert(0);
  return -1;
}

/* 
   this generic function retrieves the data type
   from the file and uses it to call a specialized 
   function.
*/

int taucs_ooc_solve_llt (void* L /* actual type: taucs_io_handle* */,
			 void* x, void* b)
{
  int flags;

  taucs_io_read((taucs_io_handle*)L,
		6,1,1,TAUCS_INT,
		&flags);

#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    return taucs_dooc_solve_llt(L,x,b);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    return taucs_sooc_solve_llt(L,x,b);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zooc_solve_llt(L,x,b);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cooc_solve_llt(L,x,b);
#endif

  assert(0);
  return -1;
}

#endif /* TAUCS_CORE_GENRAL */

/*************************************************************/
/* end of file                                               */
/*************************************************************/




