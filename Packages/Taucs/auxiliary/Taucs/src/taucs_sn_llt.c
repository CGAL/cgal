/*************************************************************/
/*                                                           */
/*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NDEBUG
#include <assert.h>

#define TAUCS_CORE_CILK
#include "taucs.h"

#ifdef TAUCS_CILK
#pragma lang -C
#endif

#ifndef TAUCS_CORE_GENERAL
#ifdef TAUCS_CILK

/*** GEMM ***/

#define TAUCS_THRESHOLD_GEMM_SMALL 20
#define TAUCS_THRESHOLD_GEMM_BLAS  80

static void
taucs_gemm_NCm1p1_small(int m, int n, int k, 
			taucs_datatype* A, int lda,
			taucs_datatype* B, int ldb,
			taucs_datatype* C, int ldc)
{
  int j,i,l;
  taucs_datatype* Cj;
  taucs_datatype* Ail;
  taucs_datatype* Bjl;
  taucs_datatype  Cij;

  Cj = C;
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      Cij = *Cj;
      Ail = A + i;
      Bjl = B + j;
      for (l=0; l<k; l++) {
	Cij = taucs_sub( Cij, taucs_mul( *Ail, taucs_conj( *Bjl ) ) );
	Ail += lda;
	Bjl += ldb;
      }
      *Cj = Cij;
      Cj++;
    }
    Cj = Cj + ldc - m;
    /* now Cj is at the top of column j+1 */
  }
}

cilk static void 
taucs_cilk_gemm(char* transa, char* transb,
		int* pm, int*  pn, int* pk,
		taucs_real_datatype* alpha, 
		taucs_datatype *A, int *plda,
		taucs_datatype *B, int *pldb,
		taucs_real_datatype* beta, 
		taucs_datatype *C, int *pldc)
{
  int    m  = *pm;
  int    n  = *pn;
  int    k  = *pk;

  assert(*transa == 'N');
  assert(*transb == 'C');
  assert(*alpha  ==-1.0);
  assert(*beta   == 1.0);

  if (n <= TAUCS_THRESHOLD_GEMM_SMALL && k <= TAUCS_THRESHOLD_GEMM_SMALL) {
    /*fprintf(stderr,"GEMM SMALL\n");*/
    taucs_gemm_NCm1p1_small(m,n,k,A,*plda,B,*pldb,C,*pldc);
    return;
  }

  if (n <= TAUCS_THRESHOLD_GEMM_BLAS && k <= TAUCS_THRESHOLD_GEMM_BLAS) {
    /*fprintf(stderr,"GEMM BLAS\n");*/
    taucs_gemm(transa, transb,
	       pm, pn, pk,
	       alpha, 
	       A, plda,
	       B, pldb,
	       beta,
	       C, pldc);
    return;
  }

  if (k >= n && k >= m) {
    int khalf1 = k/2;
    int khalf2 = k-khalf1;
    int lda = *plda;
    int ldb = *pldb;
    /*fprintf(stderr,"GEMM K/2\n");*/
    spawn taucs_cilk_gemm(transa,transb, 
			  pm, pn, &khalf1,
			  alpha, 
			  A, plda,
			  B, pldb,
			  beta,
			  C, pldc);
    sync;
    spawn taucs_cilk_gemm(transa,transb, 
			  pm, pn, &khalf2,
			  alpha, 
			  A + khalf1*lda, plda,
			  B + khalf1*ldb, pldb,
			  beta,
			  C, pldc);
    sync;
    return;
  } 

  if (n >= k && n >= m) {
    int nhalf1 = n/2;
    int nhalf2 = n-nhalf1;
    int ldc = *pldc;
    /*fprintf(stderr,"GEMM N/2\n");*/

    spawn taucs_cilk_gemm(transa,transb, 
			  pm, &nhalf1, pk,
			  alpha, 
			  A, plda,
			  B, pldb,
			  beta,
			  C, pldc);


    spawn taucs_cilk_gemm(transa,transb, 
			  pm, &nhalf2, pk,
			  alpha, 
			  A, plda,
			  B + nhalf1, pldb,
			  beta,
			  C + nhalf1*ldc, pldc);
    sync;
    return;
  }

  if (1 /* m >= k && m >= n*/) { /* the condition must be true */
    int mhalf1 = m/2;
    int mhalf2 = m-mhalf1;
    /*fprintf(stderr,"GEMM M/2\n");*/

    spawn taucs_cilk_gemm(transa,transb, 
			  &mhalf1, pn, pk,
			  alpha, 
			  A, plda,
			  B, pldb,
			  beta,
			  C, pldc);

    spawn taucs_cilk_gemm(transa,transb, 
			  &mhalf2, pn, pk,
			  alpha, 
			  A + mhalf1, plda,
			  B, pldb,
			  beta,
			  C + mhalf1, pldc);
    sync;
    return;
  }

  assert(0);

}

/*** HERK ***/

#define TAUCS_THRESHOLD_HERK_SMALL 20
#define TAUCS_THRESHOLD_HERK_BLAS  80

static void
taucs_herk_LNm1p1_small(int n, int k, 
			taucs_datatype* A, int lda,
			taucs_datatype* C, int ldc)
{
  int j,i,l;
  taucs_datatype* Cj;
  taucs_datatype* Ail;
  taucs_datatype* Ajl;
  taucs_datatype  Cij;

  Cj = C;
  for (j=0; j<n; j++) {
    for (i=j; i<n; i++) {
      Cij = *Cj;
      Ail = A + i;
      Ajl = A + j;
      for (l=0; l<k; l++) {
	Cij = taucs_sub( Cij, taucs_mul( *Ail, taucs_conj( *Ajl ) ) );
	Ail += lda;
	Ajl += lda;
      }
      *Cj = Cij;
      Cj++;
    }
    Cj = Cj + ldc - n;
    /* now Cj is at the top of column j+1, move to the diagonal */
    Cj = Cj + j+1;
  }
}

cilk static void 
taucs_cilk_herk(char* uplo, char* trans,
		int*  pn, int* pk,
		taucs_real_datatype* alpha, 
		taucs_datatype *A, int *plda,
		taucs_real_datatype* beta, 
		taucs_datatype *C, int *pldc)
{
  int    n  = *pn;
  int    k  = *pk;

  assert(*uplo  == 'L');
  assert(*trans == 'N');
  assert(*alpha ==-1.0);
  assert(*beta  == 1.0);

  if (n <= TAUCS_THRESHOLD_HERK_SMALL && k <= TAUCS_THRESHOLD_HERK_SMALL) {
    /*fprintf(stderr,"HERK SMALL\n");*/
    taucs_herk_LNm1p1_small(n,k,A,*plda,C,*pldc);
    return;
  }

  if (n <= TAUCS_THRESHOLD_HERK_BLAS && k <= TAUCS_THRESHOLD_HERK_BLAS) {
    /*fprintf(stderr,"HERK BLAS\n");*/
    taucs_herk(uplo,trans,
	       pn, pk,
	       alpha, 
	       A, plda,
	       beta,
	       C, pldc);
    return;
  }

  if (k > n) {
    int khalf1 = k/2;
    int khalf2 = k-khalf1;
    int lda = *plda;
    /*fprintf(stderr,"HERK K/2\n");*/
    spawn taucs_cilk_herk(uplo,trans, 
			  pn, &khalf1,
			  alpha, 
			  A, plda,
			  beta,
			  C, pldc);
    sync;
    spawn taucs_cilk_herk(uplo,trans,
			  pn, &khalf2,
			  alpha, 
			  A + khalf1*lda, plda,
			  beta,
			  C, pldc);
    sync;
    return;
  } else {
    int ldc = *pldc;
    int nhalf1 = n/2;
    int nhalf2 = n-nhalf1;
    /*fprintf(stderr,"HERK N/2\n");*/
    spawn taucs_cilk_herk(uplo,trans, 
			  &nhalf1, pk,
			  alpha, 
			  A, plda,
			  beta,
			  C, pldc);

    spawn taucs_cilk_gemm("No Transpose", "Conjugate", 
			  &nhalf2, &nhalf1, pk, 
			  &taucs_minusone_const, 
			  A+nhalf1  , plda,
			  A         , plda,
			  &taucs_one_const, 
			  C +nhalf1, pldc);

    spawn taucs_cilk_herk(uplo,trans, 
			  &nhalf2, pk,
			  alpha, 
			  A + nhalf1, plda,
			  beta,
			  C + nhalf1*ldc + nhalf1, pldc);

    sync;
    return;
  }

}

/*** TRSM ***/

#define TAUCS_THRESHOLD_TRSM_SMALL 20
#define TAUCS_THRESHOLD_TRSM_BLAS  80

static void
taucs_trsm_RLCNp1_small(int m, int n, 
			taucs_datatype* A, int lda,
			taucs_datatype* B, int ldb)
{
  int j,i,k;
  taucs_datatype* Bi;
  taucs_datatype* Bik;
  taucs_datatype* Ajk;
  taucs_datatype  Bij;

  Bi = B;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      Bij = *(Bi + j*ldb);
      Bik = Bi;
      Ajk = A + j;
      for (k=0; k<j; k++) {
	Bij = taucs_sub( Bij, taucs_mul( *Bik, taucs_conj( *Ajk ) ) );
	Bik += ldb;
	Ajk += lda;
      }
      *(Bi + j*ldb) = taucs_div( Bij, taucs_conj( *Ajk ) ); 
    }
    Bi = Bi + 1;
  }
}

cilk static void 
taucs_cilk_trsm(char* side, char* uplo, char* transa, char* diag, 
		int*  pm, int* pn,
		taucs_datatype* alpha, 
		taucs_datatype *A, int *plda,
		taucs_datatype *B, int *pldb)
{
  int    n  = *pn;
  int    m  = *pm;

  assert(*side   == 'R');
  assert(*uplo   == 'L');
  assert(*transa == 'C');
  assert(*diag   == 'N');
  assert(taucs_re(*alpha) == 1.0);
  assert(taucs_im(*alpha) == 0.0);

  if (m <= TAUCS_THRESHOLD_TRSM_SMALL && n <= TAUCS_THRESHOLD_TRSM_SMALL) {
    /*fprintf(stderr,"TRSM SMALL\n");*/
    taucs_trsm_RLCNp1_small(m,n,A,*plda,B,*pldb);
    return;
  }

  if (m <= TAUCS_THRESHOLD_TRSM_BLAS && n <= TAUCS_THRESHOLD_TRSM_BLAS) {
    /*fprintf(stderr,"TRSM BLAS\n");*/
    taucs_trsm(side,uplo,transa,diag, 
	       pm, pn,
	       alpha, 
	       A, plda,
	       B, pldb);
    return;
  }

  if (m >= n) {
    int mhalf1 = m/2;
    int mhalf2 = m-mhalf1;
    /*fprintf(stderr,"TRSM M/2\n");*/
    spawn taucs_cilk_trsm(side,uplo,transa,diag, 
			  &mhalf1, pn,
			  alpha, 
			  A, plda,
			  B, pldb);
    spawn taucs_cilk_trsm(side,uplo,transa,diag, 
			  &mhalf2, pn,
			  alpha, 
			  A, plda,
			  B+mhalf1, pldb);
    sync;
    return;
  } else {
    int lda = *plda;
    int ldb = *pldb;
    int nhalf1 = n/2;
    int nhalf2 = n-nhalf1;
    /*fprintf(stderr,"TRSM N/2\n");*/
    spawn taucs_cilk_trsm(side,uplo,transa,diag, 
			  pm, &nhalf1,
			  alpha, 
			  A, plda,
			  B, pldb);
    sync;
    spawn taucs_cilk_gemm("No Transpose", "Conjugate", 
			  pm, &nhalf2, &nhalf1, 
			  &taucs_minusone_const, 
			  B                  , pldb,
			  A+nhalf1           , plda,
			  &taucs_one_const, 
			  B       +nhalf1*ldb, pldb);
    sync;
    spawn taucs_cilk_trsm(side,uplo,transa,diag, 
			  pm, &nhalf2,
			  alpha, 
			  A+nhalf1+nhalf1*lda, plda,
			  B       +nhalf1*ldb, pldb);
    sync;
    return;
  }
}

/*** POTRF ***/

#define TAUCS_THRESHOLD_POTRF_SMALL 20
#define TAUCS_THRESHOLD_POTRF_BLAS  80

/*
  this routine is for Lower only, and returns
  0 or the index of the nonpositive diagonal.
*/
static int
taucs_potrf_lower_small(int n, taucs_datatype* A, int lda)
{
  int j,i,k;
  taucs_datatype* Aj;
  taucs_datatype* Ajk;
  taucs_datatype* Aik;
  taucs_datatype  Aij;
  taucs_datatype  scale;

  Aj = A;
  for (j=0; j<n; j++) {
    for (i=j; i<n; i++) {
      Aij = taucs_zero_const;
      Ajk = A+j; /* k = 0 */
      Aik = A+i; /* k = 0 */
      for (k=0; k<j; k++) {
	Aij = taucs_add( Aij , taucs_mul( (*Aik) , (*Ajk) ) );
	Aik += lda;
	Ajk += lda;
      }
      Aj[i] = taucs_sub( Aj[i] , Aij );
    }
    if ( taucs_re(Aj[j]) < 0 ) return j+1;
    scale = taucs_div( taucs_one_const, taucs_sqrt(Aj[j]) );
    for (i=j; i<n; i++) 
      Aj[i] = taucs_mul(Aj[i] , scale);
    Aj += lda;
  }

  return 0;
}

cilk static void 
taucs_cilk_potrf(char* uplo, 
		 int*  pn,
		 taucs_datatype* A, int* plda,
		 int*  pinfo)
{
  int    n  = *pn;
  int nhalf1,nhalf2;

  assert(*uplo == 'L');

  if (n <= TAUCS_THRESHOLD_POTRF_SMALL) {
    /*fprintf(stderr,"POTRF SMALL\n");*/
    *pinfo = taucs_potrf_lower_small(*pn,A,*plda);
    return;
  }

  if (n <= TAUCS_THRESHOLD_POTRF_BLAS) {
    /*fprintf(stderr,"POTRF BLAS\n");*/
    taucs_potrf(uplo,pn,A,plda,pinfo);
    return;
  }

  /*fprintf(stderr,"POTRF RECIRSIVE\n");*/

  nhalf1 = n/2;
  nhalf2 = n-nhalf1;

  spawn taucs_cilk_potrf(uplo,&nhalf1,A,plda,pinfo);
  sync;

  if (*pinfo) return;

  spawn taucs_cilk_trsm ("Right",
			 "Lower",
			 "Conjugate",
			 "No unit diagonal",
			 &nhalf2,&nhalf1,
			 &taucs_one_const,
			 A,       plda,
			 A+nhalf1,plda);
  sync;
  /*  
  taucs_trsm ("Right",
	      "Lower",
	      "Conjugate",
	      "No unit diagonal",
	      &nhalf2,&nhalf1,
	      &taucs_one_const,
	      A,       plda,
	      A+nhalf1,plda);
  */

  spawn taucs_cilk_herk ("Lower",
			 "No Conjugate",
			 &nhalf2,&nhalf1,
			 &taucs_minusone_real_const,
			 A+nhalf1,plda,
			 &taucs_one_real_const,
			 A+nhalf1+(nhalf1 * *plda), plda);
  sync;

  spawn taucs_cilk_potrf(uplo,&nhalf2,A+nhalf1+(nhalf1 * *plda),plda,pinfo);
  sync;

  if (*pinfo) *pinfo += nhalf1;
}
#else
#define taucs_cilk_potrf taucs_potrf
#define taucs_cilk_gemm  taucs_gemm
#define taucs_cilk_trsm  taucs_trsm
#define taucs_cilk_herk  taucs_herk
#endif
#endif

/*************************************************************/
/* These are really generic routines                         */
/*************************************************************/

#if 0
#ifdef TAUCS_CORE_GENERAL
void* taucs_cilk_init() {
#ifdef TAUCS_CILK
  CilkContext* context;
  int argc;
#define CILK_ACTIVE_SIZE 8
  char* argv[] = {"program_name","--nproc","8",0};
  
  for (argc=0; argv[argc]; argc++);

  taucs_printf("taucs_cilk_init\n");
  context = Cilk_init(&argc,argv);
  return context;
#else
  taucs_printf("taucs_cilk_init: This is not a Cilk build\n");
  return NULL;
#endif
}

void taucs_cilk_terminate(void* context) {
#ifdef TAUCS_CILK
  Cilk_terminate((CilkContext*) context);
#endif
}

#endif /* TAUCS_CORE_GENERAL */
#endif
/*************************************************************/
/* End of Cilk-related generic routines                      */
/*************************************************************/

#if 0
/*omer added this, I don't know why yet*/
#ifdef TAUCS_CORE_GENERAL
taucs_double taucs_dzero_const     =  0.0;
taucs_double taucs_done_const      =  1.0;
taucs_double taucs_dminusone_const = -1.0;

taucs_single taucs_szero_const     =  0.0f;
taucs_single taucs_sone_const      =  1.0f;
taucs_single taucs_sminusone_const = -1.0f;
#endif
#endif 

/*************************************************************/
/* structures                                                */
/*************************************************************/

#define FALSE 0
#define TRUE  1

/*#define BLAS_FLOPS_CUTOFF  1000.0*/
#define BLAS_FLOPS_CUTOFF  -1.0
#define SOLVE_DENSE_CUTOFF 5

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

#define SFM_F1 f1
#define SFM_F2 f2
#define SFM_U   u

typedef struct {
  int     flags;

  char    uplo;     /* 'u' for upper, 'l' for lower, ' ' don't know; prefer lower. */
  int     n;        /* size of matrix */
  int     n_sn;     /* number of supernodes */

  int* parent;      /* supernodal elimination tree */
  int* first_child; 
  int* next_child;

  int* sn_size;     /* size of supernodes (diagonal block) */
  int* sn_up_size;  /* size of subdiagonal update blocks   */
  int** sn_struct;  /* row structure of supernodes         */

  int* sn_blocks_ld;  /* lda of supernode blocks */
  taucs_datatype** sn_blocks; /* supernode blocks        */
    
  int* up_blocks_ld;  /* lda of update blocks    */
  taucs_datatype** up_blocks; /* update blocks           */
} supernodal_factor_matrix;

#ifdef TAUCS_CORE_GENERAL
/*************************************************************/
/* for qsort                                                 */
/*************************************************************/

/* this is never used */
/*
static int compare_ints(void* vx, void* vy)
{
  int* ix = (int*)vx;
  int* iy = (int*)vy;
  if (*ix < *iy) return -1;
  if (*ix > *iy) return  1;
  return 0;
}
*/

static int* compare_indirect_map;
static int compare_indirect_ints( const void* vx, const void* vy)
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

#endif /* TAUCS_CORE_GENERAL */
/*************************************************************/
/* create and free the factor object                         */
/*************************************************************/

#ifndef TAUCS_CORE_GENERAL

static supernodal_factor_matrix*
multifrontal_supernodal_create()
{
  supernodal_factor_matrix* L;
  
  L = (supernodal_factor_matrix*) taucs_malloc(sizeof(supernodal_factor_matrix));
  if (!L) return NULL;

#ifdef TAUCS_CORE_SINGLE
  L->flags = TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DOUBLE
  L->flags = TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
  L->flags = TAUCS_SCOMPLEX;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  L->flags = TAUCS_DCOMPLEX;
#endif

  L->uplo      = 'l';
  L->n         = -1; /* unused */

  L->sn_struct   = NULL;
  L->sn_size     = NULL;
  L->sn_up_size  = NULL;
  L->parent      = NULL;
  L->first_child = NULL;
  L->next_child  = NULL;
  L->sn_blocks_ld  = NULL;
  L->sn_blocks     = NULL;
  L->up_blocks_ld  = NULL;
  L->up_blocks     = NULL;

  return L;
}

void taucs_dtl(supernodal_factor_free)(void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int sn;

  if (!L) return;
  
  taucs_free(L->parent);
  taucs_free(L->first_child);
  taucs_free(L->next_child);

  taucs_free(L->sn_size);
  taucs_free(L->sn_up_size);
  taucs_free(L->sn_blocks_ld);
  taucs_free(L->up_blocks_ld);

  if (L->sn_struct)   
    for (sn=0; sn<L->n_sn; sn++)
      taucs_free(L->sn_struct[sn]);

  if (L->sn_blocks)   
    for (sn=0; sn<L->n_sn; sn++)
      taucs_free(L->sn_blocks[sn]);

  if (L->up_blocks)   
    for (sn=0; sn<L->n_sn; sn++)
      taucs_free(L->up_blocks[sn]);

  taucs_free(L->sn_struct);
  taucs_free(L->sn_blocks);
  taucs_free(L->up_blocks);

  taucs_free(L);
}

void taucs_dtl(supernodal_factor_free_numeric)(void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int sn;
  
  for (sn=0; sn<L->n_sn; sn++) {
    taucs_free(L->sn_blocks[sn]);
    L->sn_blocks[sn] = NULL;
    taucs_free(L->up_blocks[sn]);
    L->up_blocks[sn] = NULL;
  }
}

taucs_ccs_matrix*
taucs_dtl(supernodal_factor_to_ccs)(void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  taucs_ccs_matrix* C;
  int n,nnz;
  int i,j,ip,jp,sn,next;
  taucs_datatype v;
  int* len;

  n = L->n;

  len = (int*) taucs_malloc(n*sizeof(int));
  if (!len) return NULL;

  nnz = 0;
  /*
  for (sn=0; sn<L->n_sn; sn++) {
    for (jp=0; jp<(L->sn_size)[sn]; jp++) {
      j = (L->sn_struct)[sn][jp];
      len[j] = (L->sn_up_size)[sn] - jp;
      nnz += len[j];
    }
  }
  */

  for (sn=0; sn<L->n_sn; sn++) {
    for (jp=0; jp<(L->sn_size)[sn]; jp++) {
      j = (L->sn_struct)[sn][jp];
      len[j] = 0;

      for (ip=jp; ip<(L->sn_size)[sn]; ip++) {
	i = (L->sn_struct)[sn][ ip ];
	v = (L->sn_blocks)[sn][ jp*(L->sn_blocks_ld)[sn] + ip ];

	if (taucs_re(v) || taucs_im(v)) { 
	  len[j] ++;
	  nnz ++;
	}
      }
      for (ip=(L->sn_size)[sn]; ip<(L->sn_up_size)[sn]; ip++) {
	i = (L->sn_struct)[sn][ ip ];
	v = (L->up_blocks)[sn][ jp*(L->up_blocks_ld)[sn] + (ip-(L->sn_size)[sn]) ];

	if (taucs_re(v) || taucs_im(v)) { 
	  len[j] ++;
	  nnz ++;
	}
      }
    }
  }


  C = taucs_dtl(ccs_create)(n,n,nnz);
  if (!C) {
    taucs_free(len);
    return NULL;
  }

#ifdef TAUCS_CORE_SINGLE
  C->flags = TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DOUBLE
  C->flags = TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
  C->flags = TAUCS_SCOMPLEX;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  C->flags = TAUCS_DCOMPLEX;
#endif

  C->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

  (C->colptr)[0] = 0;
  for (j=1; j<=n; j++) (C->colptr)[j] = (C->colptr)[j-1] + len[j-1];

  taucs_free(len);

  for (sn=0; sn<L->n_sn; sn++) {
    for (jp=0; jp<(L->sn_size)[sn]; jp++) {
      j = (L->sn_struct)[sn][jp];

      next = (C->colptr)[j];

      /*
      memcpy((C->rowind) + next,
	     ((L->sn_struct)[sn]) + jp,
	     ((L->sn_up_size)[sn] - jp) * sizeof(int));
      memcpy((C->taucs_values) + next,
	     ((L->sn_blocks)[sn]) + (jp*(L->sn_blocks_ld)[sn] + jp),
	     ((L->sn_size)[sn] - jp) * sizeof(taucs_datatype));
      next += ((L->sn_size)[sn] - jp);
      memcpy((C->taucs_values) + next,
	     ((L->up_blocks)[sn]) + jp*(L->up_blocks_ld)[sn],
	     ((L->sn_up_size)[sn] - (L->sn_size)[sn]) * sizeof(taucs_datatype));
      */

      for (ip=jp; ip<(L->sn_size)[sn]; ip++) {
	i = (L->sn_struct)[sn][ ip ];
	v = (L->sn_blocks)[sn][ jp*(L->sn_blocks_ld)[sn] + ip ];

	if (!taucs_re(v) && !taucs_im(v)) continue;
	/*if (v == 0.0) continue;*/

	(C->rowind)[next] = i;
	(C->taucs_values)[next] = v;
	next++;
      }
      for (ip=(L->sn_size)[sn]; ip<(L->sn_up_size)[sn]; ip++) {
	i = (L->sn_struct)[sn][ ip ];
	v = (L->up_blocks)[sn][ jp*(L->up_blocks_ld)[sn] + (ip-(L->sn_size)[sn]) ];

	if (!taucs_re(v) && !taucs_im(v)) continue;
	/*if (v == 0.0) continue;*/

	(C->rowind)[next] = i;
	(C->taucs_values)[next] = v;
	next++;
      }
    }
  }

  return C;
}

/* just get the diagonal of a supernodal factor, for Penny */

taucs_datatype*
taucs_dtl(supernodal_factor_get_diag)(void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int j,ip,jp,sn;/*i,next omer*/
  taucs_datatype  v;
  taucs_datatype* diag;

  diag = (taucs_datatype*) taucs_malloc((L->n) * sizeof(taucs_datatype));
  if (!diag) return NULL;

  for (sn=0; sn<L->n_sn; sn++) {
    for (jp=0; jp<(L->sn_size)[sn]; jp++) {
      j = (L->sn_struct)[sn][jp];

      ip=jp; /* we just want the diagonal */
      
      v = (L->sn_blocks)[sn][ jp*(L->sn_blocks_ld)[sn] + ip ];
      
      diag[ j ] = v;
    }
  }

  return diag;
}


/*************************************************************/
/* create and free frontal matrices                          */
/*************************************************************/

static supernodal_frontal_matrix* 
supernodal_frontal_create(int* firstcol_in_supernode,
			  int sn_size,
			  int n, 
			  int* rowind)
{
  supernodal_frontal_matrix* tmp;

  tmp = (supernodal_frontal_matrix*)taucs_malloc(sizeof(supernodal_frontal_matrix));
  if(tmp==NULL) return NULL;

  tmp->sn_size = sn_size;
  tmp->n = n;

  tmp->rowind = rowind;

  tmp->n = n;
  tmp->sn_size = sn_size;
  tmp->up_size = n-sn_size;

  tmp->sn_vertices = rowind;
  tmp->up_vertices = rowind + sn_size;

  /* on some platforms, malloc(0) fails, so we avoid such calls */

  tmp->SFM_F1 = tmp->SFM_F2 = tmp->SFM_U = NULL;

  if (tmp->sn_size)
    tmp->SFM_F1 = (taucs_datatype*)taucs_calloc((tmp->sn_size)*(tmp->sn_size),sizeof(taucs_datatype));

  if (tmp->sn_size && tmp->up_size)
    tmp->SFM_F2 = (taucs_datatype*)taucs_calloc((tmp->up_size)*(tmp->sn_size),sizeof(taucs_datatype));

  if (tmp->up_size)
    tmp->SFM_U  = (taucs_datatype*)taucs_calloc((tmp->up_size)*(tmp->up_size),sizeof(taucs_datatype));

  if((   tmp->SFM_F1==NULL && tmp->sn_size)
     || (tmp->SFM_F2==NULL && tmp->sn_size && tmp->up_size)
     || (tmp->SFM_U ==NULL && tmp->up_size)) {
    taucs_free(tmp->SFM_U);
    taucs_free(tmp->SFM_F1);
    taucs_free(tmp->SFM_F2);
    taucs_free(tmp);
    return NULL;
  }

  assert(tmp);
  return tmp;
}

static void supernodal_frontal_free(supernodal_frontal_matrix* to_del)
{
  /* 
     SFM_F1 and SFM_F2 are moved to the factor,
     but this function may be called before they are
     moved.
  */


  if (to_del) {
    taucs_free(to_del->SFM_F1);
    taucs_free(to_del->SFM_F2);
    taucs_free(to_del->SFM_U);
    taucs_free(to_del);
  }
}

/*************************************************************/
/* factor a frontal matrix                                   */
/*************************************************************/

cilk
static int
multifrontal_supernodal_front_factor(int sn,
				     int* firstcol_in_supernode,
				     int sn_size,
				     taucs_ccs_matrix* A,
				     supernodal_frontal_matrix* mtr,
				     int* bitmap,
				     supernodal_factor_matrix* snL)
{
  int i,j;
  int* ind;
  taucs_datatype* re;
  int INFO;

  /* creating transform for real indices */
  for(i=0;i<mtr->sn_size;i++) bitmap[mtr->sn_vertices[i]] = i;
  for(i=0;i<mtr->up_size;i++) bitmap[mtr->up_vertices[i]] = mtr->sn_size + i;

  /* adding sn_size column of A to first sn_size column of frontal matrix */

  for(j=0;j<(mtr->sn_size);j++) {
    ind = &(A->rowind[A->colptr[*(firstcol_in_supernode+j)]]);
    re  = &(A->taucs_values[A->colptr[*(firstcol_in_supernode+j)]]); 
    for(i=0;
	i < A->colptr[*(firstcol_in_supernode+j)+1] 
            - A->colptr[*(firstcol_in_supernode+j)];
	i++) {
      if (bitmap[ind[i]] < mtr->sn_size)
	mtr->SFM_F1[ (mtr->sn_size)*j + bitmap[ind[i]]] =
	  taucs_add( mtr->SFM_F1[ (mtr->sn_size)*j + bitmap[ind[i]]] , re[i] );
      else
	mtr->SFM_F2[ (mtr->up_size)*j + bitmap[ind[i]] - mtr->sn_size] =
	  taucs_add( mtr->SFM_F2[ (mtr->up_size)*j + bitmap[ind[i]] - mtr->sn_size] , re[i] );
    }
  }

  /* we use the BLAS through the Fortran interface */

  /* solving of lower triangular system for L */
  if (mtr->sn_size) {
    /*
    taucs_potrf ("LOWER",
		 &(mtr->sn_size),
		 mtr->SFM_F1,&(mtr->sn_size),
		 &INFO);
    */
    spawn taucs_cilk_potrf ("LOWER",
		 &(mtr->sn_size),
		 mtr->SFM_F1,&(mtr->sn_size),
		 &INFO);
    sync;
  }


  if (INFO) {
    taucs_printf("sivan %d %d\n",sn,sn_size);
    taucs_printf("\t\tLL^T Factorization: Matrix is not positive definite.\n");
    taucs_printf("\t\t                    nonpositive pivot in column %d\n",
		 mtr->sn_vertices[INFO-1]);
    return -1;
  }

  /* getting completion for found columns of L */
  if (mtr->up_size && mtr->sn_size) {

    spawn taucs_cilk_trsm ("Right",
			   "Lower",
			   "Conjugate",
			   "No unit diagonal",
			   &(mtr->up_size),&(mtr->sn_size),
			   &taucs_one_const,
			   mtr->SFM_F1,&(mtr->sn_size),
			   mtr->SFM_F2,&(mtr->up_size));
    sync;
    /*
    taucs_trsm ("Right",
		"Lower",
		"Conjugate",
		"No unit diagonal",
		&(mtr->up_size),&(mtr->sn_size),
		&taucs_one_const,
		mtr->SFM_F1,&(mtr->sn_size),
		mtr->SFM_F2,&(mtr->up_size));
    */
  }

  (snL->sn_blocks   )[sn] = mtr->SFM_F1;
  (snL->sn_blocks_ld)[sn] = mtr->sn_size;

  (snL->up_blocks   )[sn] = mtr->SFM_F2;
  (snL->up_blocks_ld)[sn] = mtr->up_size;
  /* printf("*** sn=%d up_ld=%d (%d)\n",sn,mtr->up_size,(snL->up_vertex_ptr)[sn+1] - (snL->up_vertex_ptr)[sn]);*/

  /* computation of updated part of frontal matrix */
  if (mtr->up_size && mtr->sn_size) {
    spawn taucs_cilk_herk ("Lower",
			   "No Conjugate",
			   &(mtr->up_size),&(mtr->sn_size),
			   &taucs_minusone_real_const,
			   mtr->SFM_F2,&(mtr->up_size),
			   &taucs_one_real_const,
			   mtr->SFM_U, &(mtr->up_size));
    sync;
  }

  mtr->SFM_F1 = NULL; /* so we don't free twice */
  mtr->SFM_F2 = NULL; /* so we don't free twice */

  return 0;
 }

/*************************************************************/
/* extend-add                                                */
/*************************************************************/

static void 
multifrontal_supernodal_front_extend_add(
					 supernodal_frontal_matrix* parent_mtr,
					 supernodal_frontal_matrix* my_mtr,
					 int* bitmap)
{
  int j,i,parent_i,parent_j;
  taucs_datatype v;

  for(i=0;i<parent_mtr->sn_size;i++) bitmap[parent_mtr->sn_vertices[i]] = i;
  for(i=0;i<parent_mtr->up_size;i++) bitmap[parent_mtr->up_vertices[i]] = (parent_mtr->sn_size)+i;

  /* extend add operation for update matrix */
  for(j=0;j<my_mtr->up_size;j++) {
    for(i=j;i<my_mtr->up_size;i++) {
      parent_j = bitmap[ my_mtr->up_vertices[j] ];
      parent_i = bitmap[ my_mtr->up_vertices[i] ];
      /* we could skip this if indices were sorted */
      if (parent_j>parent_i) {
	int tmp = parent_j;
	parent_j = parent_i;
	parent_i = tmp;
      }

      v = (my_mtr->SFM_U)[(my_mtr->up_size)*j+i];

      if (parent_j < parent_mtr->sn_size) {
	if (parent_i < parent_mtr->sn_size) {
	  (parent_mtr->SFM_F1)[ (parent_mtr->sn_size)*parent_j + parent_i] =
	    taucs_add( (parent_mtr->SFM_F1)[ (parent_mtr->sn_size)*parent_j + parent_i] , v );
	} else {
	  (parent_mtr->SFM_F2)[ (parent_mtr->up_size)*parent_j + (parent_i-parent_mtr->sn_size)] =
	    taucs_add( (parent_mtr->SFM_F2)[ (parent_mtr->up_size)*parent_j + (parent_i-parent_mtr->sn_size)] , v );
	}
      } else {
	(parent_mtr->SFM_U)[ (parent_mtr->up_size)*(parent_j-parent_mtr->sn_size) + (parent_i-parent_mtr->sn_size)] =
	  taucs_add( (parent_mtr->SFM_U)[ (parent_mtr->up_size)*(parent_j-parent_mtr->sn_size) + (parent_i-parent_mtr->sn_size)] , v);
      }
    }
  }
}

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*************************************************************/
/* symbolic elimination                                      */
/*************************************************************/

#ifdef TAUCS_CORE_GENERAL

/* UNION FIND ROUTINES */

static int uf_makeset(int* uf, int i)        { uf[i] = i; return i; }
static int uf_find   (int* uf, int i)
{ 
  if (uf[i] != i) 
    uf[i] = uf_find(uf,uf[i]); 
  return uf[i]; 
}
static int uf_union  (int* uf, int s, int t) {
  if (uf_find(uf,s) < uf_find(uf,t)) {
    uf[uf_find(uf,s)] = uf_find(uf,t); 
    return (uf_find(uf,t)); 
  } else {
    uf[uf_find(uf,s)] = uf_find(uf,t); 
    return (uf_find(uf,t)); 
  }
}

static
void recursive_postorder(int  j,
			 int  first_child[],
			 int  next_child[],
			 int  postorder[],
			 int  ipostorder[],
			 int* next)
{
  int c;
  for (c=first_child[j]; c != -1; c = next_child[c]) {
    /*printf("*** %d is child of %d\n",c,j);*/
    recursive_postorder(c,first_child,next_child,
			postorder,ipostorder,next);
  }
  /*printf(">>> j=%d next=%d\n",j,*next);*/
  if (postorder)  postorder [*next] = j;
  if (ipostorder) ipostorder[j] = *next;
  (*next)++;
}

#define GILBERT_NG_PEYTON_ANALYSIS_SUP

/* in a few tests the supernodal version seemed slower */
#undef GILBERT_NG_PEYTON_ANALYSIS_SUP

static int ordered_uf_makeset(int* uf, int i)
{ 
  uf[i] = i; 
  return i; 
}
static int ordered_uf_find   (int* uf, int i) 
{ 
  if (uf[i] != i) 
    uf[i] = uf_find(uf,uf[i]); 
  return uf[i]; 
}
static int ordered_uf_union  (int* uf, int s, int t) 
{
  assert(uf[t] == t);
  assert(uf[s] == s);
  assert(t > s);
  if (t > s) {
    uf[s] = t; 
    return t; 
  } else
    uf[t] = s; 
    return s; 
}

static void 
tree_level(int j,
	   int isroot, 
	   int first_child[],
	   int next_child[],
	   int level[],
	   int level_j)
{
  int c;
  if (!isroot) level[j] = level_j;
  for (c=first_child[j]; c != -1; c = next_child[c]) {
    tree_level(c,
	       FALSE,
	       first_child,
	       next_child,
	       level,
	       level_j+1);
  }
}

static void
tree_first_descendant(int j,
		      int isroot, 
		      int first_child[],
		      int next_child[],
		      int ipostorder[],
		      int first_descendant[])
{
  int c;
  int fd = ipostorder[j];
  for (c=first_child[j]; c != -1; c = next_child[c]) {
    tree_first_descendant(c,
			  FALSE,
			  first_child,
			  next_child,
			  ipostorder,
			  first_descendant);
    if (first_descendant[c] < fd) fd = first_descendant[c]; 
  }
  if (!isroot) first_descendant[j] = fd;
}


int
taucs_ccs_etree(taucs_ccs_matrix* A,
		int* parent,
		int* l_colcount,
		int* l_rowcount,
		int* l_nnz);

int 
taucs_ccs_etree_liu(taucs_ccs_matrix* A,
		    int* parent,
		    int* l_colcount,
		    int* l_rowcount,
		    int* l_nnz);



static int
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
			       int            ipostorder[]
			       )
{
  int  i,ip,c,c_sn;
  int  in_previous_sn;
  int  nnz = 0; /* just to suppress the warning */
  
  for (c=first_child[j]; c != -1; c = next_child[c]) {
    if (recursive_symbolic_elimination(c,A,
				       first_child,next_child,
				       n_sn,
				       sn_size,sn_up_size,sn_rowind,
				       sn_first_child,sn_next_child,
				       rowind, /* temporary */
				       column_to_sn_map,
				       map,
				       do_order,ipostorder
				       ) 
	== -1) return -1;
  }

  in_previous_sn = 1;
  if (j == A->n) 
    in_previous_sn = 0; /* this is not a real column */
  else if (first_child[j] == -1) 
    in_previous_sn = 0; /* this is a leaf */
  else if (next_child[first_child[j]] != -1) 
    in_previous_sn = 0; /* more than 1 child */
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

    return 0;
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
    if (!( sn_rowind [*n_sn] )) return -1;
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

  return 0;
}

/* count zeros and nonzeros in a supernode to compute the */
/* utility of merging fundamental supernodes.             */

typedef struct {
  double zeros;
  double nonzeros;
} znz;

static znz
recursive_amalgamate_supernodes(int           sn,
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
				int            ipostorder[]
				)
{
  int  i,ip,c_sn,gc_sn;
  /*int  i,ip,c,c_sn,gc_sn;*/
  int  nnz;
  int  nchildren /*, ichild*/; /* number of children, child index */
  znz* c_znz = NULL;
  znz  sn_znz, merged_znz;
  /*int zero_count = 0;*/
  int new_sn_size, new_sn_up_size;

  sn_znz.zeros    = 0.0;
  sn_znz.nonzeros = (double) (((sn_up_size[sn] - sn_size[sn]) * sn_size[sn]) 
                              + (sn_size[sn] * (sn_size[sn] + 1))/2);

  if (sn_first_child[sn] == -1) { /* leaf */
    return sn_znz;
  }

  nchildren = 0;
  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn])
    nchildren++;

  /*  c_znz = (znz*) alloca(nchildren * sizeof(znz));*/
  c_znz = (znz*) taucs_malloc(nchildren * sizeof(znz));
  assert(c_znz);

  /*printf("supernode %d out of %d\n",sn,*n_sn);*/

  /* merge the supernode with its children! */

  i = 0;
  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
    c_znz[i] = 
      recursive_amalgamate_supernodes(c_sn,
				      n_sn,
				      sn_size,sn_up_size,sn_rowind,
				      sn_first_child,sn_next_child,
				      rowind, /* temporary */
				      column_to_sn_map,
				      map,
				      do_order,ipostorder
				      );
    assert(c_znz[i].zeros + c_znz[i].nonzeros ==
	   (double) (((sn_up_size[c_sn] - sn_size[c_sn]) * sn_size[c_sn]) 
		     + (sn_size[c_sn] * (sn_size[c_sn] + 1))/2 ));
    i++;
  }

  merged_znz.nonzeros = sn_znz.nonzeros;
  merged_znz.zeros    = sn_znz.zeros;
                   
  for (i=0; i<nchildren; i++) {
    merged_znz.nonzeros += (c_znz[i]).nonzeros;
    merged_znz.zeros    += (c_znz[i]).zeros;
  }

  taucs_free(c_znz);

  /*  printf("supernode %d out of %d (continuing)\n",sn,*n_sn);*/

  /* should we merge the supernode with its children? */

  nnz = 0;
  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
    for (ip=0; ip<sn_size[c_sn]; ip++) {
      i = sn_rowind[c_sn][ip];
      assert( map[i] != sn );
      map[i] = sn;
      rowind[nnz] = i;
      nnz++;
    }
  }

  for (ip=0; ip<sn_size[sn]; ip++) {
    i = sn_rowind[sn][ip];
    assert( map[i] != sn );
    map[i] = sn;
    rowind[nnz] = i;
    nnz++;
  }

  new_sn_size = nnz;

  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
    for (ip=sn_size[c_sn]; ip<sn_up_size[c_sn]; ip++) {
      i = sn_rowind[c_sn][ip];
      if (map[i] != sn) { /* new row index */
	map[i] = sn;
	rowind[nnz] = i;
	nnz++;
      }
    }
  }

  for (ip=sn_size[sn]; ip<sn_up_size[sn]; ip++) {
    i = sn_rowind[sn][ip];
    if (map[i] != sn) { /* new row index */
      map[i] = sn;
      rowind[nnz] = i;
      nnz++;
    }
  }
  
  new_sn_up_size = nnz;

  if (do_order) {
    compare_indirect_map = ipostorder;
    qsort(rowind,nnz,sizeof(int),compare_indirect_ints);
  }

  /* determine whether we should merge the supernode and its children */

  {
    int n;
    double* zcount = NULL;

    n = 0;
    for (ip=0; ip<nnz; ip++) {
      i = rowind[ip];
      if (i >= n) n = i+1;
    }

    /*zcount = (double*) alloca(n * sizeof(double));*/
    zcount = (double*) taucs_malloc(n * sizeof(double));
    assert(zcount);
    
    for (ip=0; ip<new_sn_size; ip++) {
      i = rowind[ip]; assert(i<n);
      zcount[i] = (double) (ip+1);
    }
    for (ip=new_sn_size; ip<new_sn_up_size; ip++) {
      i = rowind[ip]; assert(i<n);
      zcount[i] = (double) new_sn_size;
    }

    /*
    for (ip=0; ip<new_sn_up_size; ip++) 
      printf("row %d zcount = %.0f\n",rowind[ip],zcount[rowind[ip]]);
    */
    
    for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
      for (ip=0; ip<sn_size[c_sn]; ip++) {
	i = sn_rowind[c_sn][ip]; assert(i<n);
	zcount[i] -= (double) (ip+1);
      }
      for (ip=sn_size[c_sn]; ip<sn_up_size[c_sn]; ip++) {
	i = sn_rowind[c_sn][ip]; assert(i<n);
	zcount[i] -= (double) sn_size[c_sn];
      }
    }

    for (ip=0; ip<sn_size[sn]; ip++) {
      i = sn_rowind[sn][ip]; assert(i<n);
      zcount[i] -= (double) (ip+1);
    }
    for (ip=sn_size[sn]; ip<sn_up_size[sn]; ip++) {
      i = sn_rowind[sn][ip]; assert(i<n);
      zcount[i] -= (double) sn_size[sn];
    }

    /*
    for (ip=0; ip<new_sn_up_size; ip++) 
      printf("ROW %d zcount = %.0f\n",rowind[ip],zcount[rowind[ip]]);
    printf("zeros before merging %.0f\n",merged_znz.zeros);
    */
    
    for (ip=0; ip<new_sn_up_size; ip++) {
      i = rowind[ip]; assert(i<n);
      assert(zcount[i] >= 0.0);
      merged_znz.zeros += zcount[i];
    }

    /*printf("zeros after merging %.0f\n",merged_znz.zeros);*/

    /* voodoo constants (need some kind of a utility function */
    if ((new_sn_size < 16)
	||
	((sn_size[sn] < 50) && (merged_znz.zeros < 0.5 * merged_znz.nonzeros))
	||
	((sn_size[sn] < 250) && (merged_znz.zeros < 0.25 * merged_znz.nonzeros))
	||
	((sn_size[sn] < 500) && (merged_znz.zeros < 0.10 * merged_znz.nonzeros))
	||
	(merged_znz.zeros < 0.05 * merged_znz.nonzeros)
	) {
      /*
      taucs_printf("merging sn %d, zeros (%f) vs nonzeros (%f)\n",
		   sn,merged_znz.zeros,merged_znz.nonzeros);
      */
    } else {
      /*
      taucs_printf("sn %d, too many zeros (%f) vs nonzeros (%f)\n",
		   sn,merged_znz.zeros,merged_znz.nonzeros);
      printf("returning without merging\n");
      */
      taucs_free(zcount);
      return sn_znz;
    }

    taucs_free(zcount);
  }

  /* now merge the children lists */

  sn_size[sn]    = new_sn_size;
  sn_up_size[sn] = new_sn_up_size;
  sn_rowind[sn]  = (int*) taucs_realloc(sn_rowind[sn], 
				  new_sn_up_size * sizeof(int));
  for (ip=0; ip<new_sn_up_size; ip++) sn_rowind[sn][ip] = rowind[ip];

  /*  printf("supernode %d out of %d (merging)\n",sn,*n_sn);*/

  nchildren = 0;
  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
    for (ip=0; ip<sn_size[c_sn]; ip++) {
      i = (sn_rowind[c_sn])[ip];
      assert(column_to_sn_map[i] == c_sn);
      column_to_sn_map[i] = sn;
    }

    for (gc_sn=sn_first_child[c_sn]; gc_sn != -1; gc_sn = sn_next_child[gc_sn]) {
      rowind[nchildren] = gc_sn;
      nchildren++;
    }
  }

  /* free the children's rowind vectors */
  for (c_sn=sn_first_child[sn]; c_sn != -1; c_sn = sn_next_child[c_sn]) {
    taucs_free( sn_rowind[c_sn] );
    sn_rowind[c_sn]  = NULL;
    sn_size[c_sn]    = 0;
    sn_up_size[c_sn] = 0;
  }

  sn_first_child[sn] = -1;
  for (i=0; i<nchildren; i++) {
    sn_next_child[ rowind[i] ] = sn_first_child[sn];
    sn_first_child[sn] = rowind[i];
  }    

  /*
  printf("supernode %d out of %d (done)\n",sn,*n_sn);
  printf("returning, merging\n");
  */
  return merged_znz;
}
#endif /* #ifdef TAUCS_CORE_GENERAL */


#ifndef TAUCS_CORE_GENERAL


/*************************************************************/
/* factor routines                                           */
/*************************************************************/

static void extend_add_wrapper(supernodal_frontal_matrix * child_matrix,
			       supernodal_frontal_matrix ** my_matrix_ptr,
			       int is_root,
			       int *v,
			       int sn_size,
			       int sn_up_size,
			       int * rowind,
			       int * bitmap,
			       int * fail) {

  if (*fail) {
    if (*my_matrix_ptr)
      supernodal_frontal_free(*my_matrix_ptr);
    return;
  }
  
  if (!is_root) {
    if (!(*my_matrix_ptr)) {
      *my_matrix_ptr = supernodal_frontal_create(v,sn_size,sn_up_size,rowind);
      if (!(*my_matrix_ptr)) {
	*fail = TRUE;
	supernodal_frontal_free(child_matrix);
	return;
      }
    }
    multifrontal_supernodal_front_extend_add(*my_matrix_ptr,child_matrix,bitmap);
  }
  
  /* moved outside "if !is_root"; Sivan 27 Feb 2002 */
  supernodal_frontal_free(child_matrix);
}

cilk
static supernodal_frontal_matrix*
recursive_multifrontal_supernodal_factor_llt(int sn,       /* this supernode */
					     int is_root,  /* is v the root? */
#if 1
					     int** bitmaps,
#else
					     int* bitmap,
#endif
					     taucs_ccs_matrix* A,
					     supernodal_factor_matrix* snL,
					     int* fail)
{
  supernodal_frontal_matrix* my_matrix=NULL;
  /*supernodal_frontal_matrix* child_matrix=NULL;*/
  int child;
  int* v; 
  int  sn_size;
  int* first_child   = snL->first_child;
  int* next_child    = snL->next_child;

#ifdef TAUCS_CILK
  /* Inlet for syncronization */
  inlet void extend_add_inlet(supernodal_frontal_matrix * child_matrix) {

    if (!is_root) {
      if (!(my_matrix)) {
	my_matrix = supernodal_frontal_create(v,sn_size,
					      snL->sn_up_size[sn],snL->sn_struct[sn]);

	if (!(my_matrix)) {
	  *fail = TRUE;
	  supernodal_frontal_free(child_matrix);
	  return;
	}
      }
      multifrontal_supernodal_front_extend_add(my_matrix,child_matrix,bitmaps[Self]);
    }

    /* moved outside "if !is_root"; Sivan 27 Feb 2002 */
    supernodal_frontal_free(child_matrix);
    /*
      The following approach is not working correctly because nothing is guaranteed about different procedure instances atomcity:
      
      extend_add_wrapper(child_matrix,&my_matrix,is_root,v,sn_size,snL->sn_up_size[sn],snL->sn_struct[sn],bitmaps[Self],fail);
    */
  }
#endif

  /* Sivan fixed a bug 25/2/2003: v was set even at the root, */
  /* but this element of sn_struct was not allocated.         */

  if (!is_root) {
    sn_size = snL->sn_size[sn];
    v = &( snL->sn_struct[sn][0] );
  } else {
    sn_size = -1;
    v = NULL; /* not used */
  }

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    /* original non-cilk code: */
    /*
    child_matrix = 
      recursive_multifrontal_supernodal_factor_llt(child,
						   FALSE,
						   bitmap,
						   A,snL,fail);
    */

#ifdef TAUCS_CILK
    extend_add_inlet(spawn recursive_multifrontal_supernodal_factor_llt(child,
									FALSE,
									bitmaps,
									A,snL,fail));
#else
    extend_add_wrapper(recursive_multifrontal_supernodal_factor_llt(child,
								    FALSE,
								    bitmaps,
								    A,snL,fail),
		       &my_matrix,
		       is_root,
		       v,
		       sn_size,
		       snL->sn_up_size[sn],
		       snL->sn_struct[sn],
		       bitmaps[Self],
		       fail);
#endif


    if (*fail) { 
      if (my_matrix) supernodal_frontal_free(my_matrix);
      return NULL;
    }

#if 0
    if (!is_root) {
      if (!my_matrix) {
	my_matrix =  supernodal_frontal_create(v,sn_size,
					       snL->sn_up_size[sn],
					       snL->sn_struct[sn]);
	if (!my_matrix) {
	  *fail = TRUE;
	  supernodal_frontal_free(child_matrix);
	  return NULL;
	}
      }
      multifrontal_supernodal_front_extend_add(my_matrix,child_matrix,bitmap);
    }
    /* moved outside "if !is_root"; Sivan 27 Feb 2002 */
    supernodal_frontal_free(child_matrix);
#endif /* 0, old pre-cilk code */
  }
  sync;

  /* in case we have no children, we allocate now */
  if (!is_root && !my_matrix) {
    my_matrix =  supernodal_frontal_create(v,sn_size,
					   snL->sn_up_size[sn],
					   snL->sn_struct[sn]);
    if (!my_matrix) {
      *fail = TRUE;
      return NULL;
    }
  }
  
  if(!is_root) {
    int rc;
    rc = spawn multifrontal_supernodal_front_factor(sn,
						    v,sn_size,
						    A,
						    my_matrix,
#if 1
						    bitmaps[Self],
#else
						    bitmap,
#endif
						    snL);
    sync;
    if (rc) { 
      /* nonpositive pivot */
      *fail = TRUE;
      supernodal_frontal_free(my_matrix);
      return NULL;
    }
  }
  return my_matrix;
}

cilk
void* 
taucs_dtl(ccs_factor_llt_mf)(taucs_ccs_matrix* A)
{
  void* p;

  p = spawn taucs_dtl(ccs_factor_llt_mf_maxdepth)(A,0);
  sync;

  return p;
}

cilk
static void
recursive_multifrontal_supernodal_factor_llt_caller(int n_sn,     /* this supernode */
						    int is_root,  /* is v the root? */
						    taucs_ccs_matrix* A,
						    supernodal_factor_matrix* snL,
						    int* fail)
{
  int** maps;
  int   i,j;
  supernodal_frontal_matrix* always_null;

  maps = (int**)taucs_malloc(Cilk_active_size*sizeof(int*));
  if (!maps) {
    taucs_supernodal_factor_free(snL);
    assert(0); return;
    /*return NULL;*/
  }

  for (i=0; i < Cilk_active_size; i++) {
    maps[i] = (int*)taucs_malloc((A->n+1)*sizeof(int));
    if (!maps[i]) {
      for (j=0; j < i ; j++)
	taucs_free(maps[j]);
      taucs_free(maps);
      taucs_supernodal_factor_free(snL);
      assert(0); return;
      /*return NULL;*/
    }
  }

  /*#ifdef TAUCS_CILK  */
#if 0
  context = Cilk_init(&argc,argv);
  always_null = EXPORT(recursive_multifrontal_supernodal_factor_llt)(context,
								     n_sn,
								     TRUE, 
								     maps,
								     A,snL,fail);
  Cilk_terminate(context);
#else
  always_null = spawn recursive_multifrontal_supernodal_factor_llt(n_sn,
								   TRUE, 
								   maps,
								   A,snL,fail);
  sync;
#endif

  for(i=0;i<Cilk_active_size;i++)
    taucs_free(maps[i]);
  taucs_free(maps);

  /*
    always_null = spawn recursive_multifrontal_supernodal_factor_llt((L->n_sn),  
    TRUE, 
    map,
    A,L,&fail);
  */
}

cilk
void* 
taucs_dtl(ccs_factor_llt_mf_maxdepth)(taucs_ccs_matrix* A,int max_depth)
{
  supernodal_factor_matrix* L;
#if 1
#else
  int* map;
#endif
  int fail;
  double wtime, ctime;

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();
  if (!L) return NULL;

#ifdef TAUCS_CORE_COMPLEX
  fail = taucs_ccs_symbolic_elimination(A,L,
					TRUE /* sort, to avoid complex conjuation */,
					max_depth);
#else
  fail = taucs_ccs_symbolic_elimination(A,L,
					FALSE /* don't sort row indices */          ,
					max_depth);
#endif
  if (fail == -1) {
    taucs_supernodal_factor_free(L);
    return NULL;
  }

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

#if 1
#else
  map = (int*)taucs_malloc((A->n+1)*sizeof(int));
  if (!map) {
    taucs_supernodal_factor_free(L);
    return NULL;
  }
#endif

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  fail = FALSE;
  spawn recursive_multifrontal_supernodal_factor_llt_caller((L->n_sn),  
							    TRUE, 
							    A,L,&fail);
  sync;
  
  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSupernodal Multifrontal LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

#if 1
#else
  taucs_free(map);
#endif

  if (fail) {
    taucs_supernodal_factor_free(L);
    return NULL;
  }

  return (void*) L;
}

/*************************************************************/
/* symbolic-numeric routines                                 */
/*************************************************************/

void* 
taucs_dtl(ccs_factor_llt_symbolic)(taucs_ccs_matrix* A)
{
  return taucs_dtl(ccs_factor_llt_symbolic_maxdepth)(A,0);
}

void* 
taucs_dtl(ccs_factor_llt_symbolic_maxdepth)(taucs_ccs_matrix* A, int max_depth)
{
  supernodal_factor_matrix* L;
  int fail;
  double wtime, ctime;

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();
  if (!L) return NULL;

#ifdef TAUCS_CORE_COMPLEX
  fail = taucs_ccs_symbolic_elimination(A,L,
					TRUE /* sort, to avoid complex conjuation */,
					max_depth);
#else
  fail = taucs_ccs_symbolic_elimination(A,L,
					FALSE /* don't sort row indices */          ,
					max_depth);
#endif

  if (fail == -1) {
    taucs_supernodal_factor_free(L);
    return NULL;
  }

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);
  return L;
}

cilk
int 
taucs_dtl(ccs_factor_llt_numeric)(taucs_ccs_matrix* A,void* vL)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  int* map;
  int fail;
  double wtime, ctime;

  map = (int*)taucs_malloc((A->n+1)*sizeof(int));

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  /* XXX: sivan, we don't need map */
  fail = FALSE;
  spawn recursive_multifrontal_supernodal_factor_llt_caller((L->n_sn),  
							    TRUE, 
							    A,L,&fail);
  sync;
  /*
    recursive_multifrontal_supernodal_factor_llt((L->n_sn),  
					       TRUE, 
					       map,
					       A,L,&fail);
  */

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSupernodal Multifrontal LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  taucs_free(map);

  if (fail) {
    taucs_supernodal_factor_free_numeric(L);
    return -1;
  }

  return 0;
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
  int first_row = 0;
  int row_count=0;
  int PK,M,N,LDA,LDB,LDC;

  /*
  for(i=0;i<sn_size_father;i++) {
    bitmap[L->sn_struct[J][i]]=i+1;
  }

  for(i=sn_size_father;i<sn_up_size_father;i++)
    bitmap[L->sn_struct[J][i]] = i - sn_size_father + 1;
  */

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
    LDA = LDB = (L->up_blocks_ld)[K];
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

#if 0
	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] -= dense_update_matrix[j*LDC+ir];
#endif

	L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] = 
	  taucs_sub( L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)] , dense_update_matrix[j*LDC+ir]);

	/*	taucs_printf("sn_block: L[%d,%d] = %lf\n",(bitmap[L->sn_struct[K][first_row+ir]]-1),(bitmap[L->sn_struct[K][first_row+j]]-1),L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*sn_size_father+(bitmap[L->sn_struct[K][first_row+ir]]-1)]);*/

      }

    for(j=0;j<row_count;j++)
      for(ir=row_count;ir<M;ir++){
#if 0
	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->up_blocks_ld)[J]+(bitmap[L->sn_struct[K][ir+first_row]]-1)] -= dense_update_matrix[j*LDC+ir];
#endif

	L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->up_blocks_ld)[J]+(bitmap[L->sn_struct[K][ir+first_row]]-1)] =
	  taucs_sub( L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->up_blocks_ld)[J]+(bitmap[L->sn_struct[K][ir+first_row]]-1)] , dense_update_matrix[j*LDC+ir]);

	/*	taucs_printf("up_block: L[%d,%d] = %lf\n",(bitmap[L->sn_struct[K][ir+first_row]]-1),(bitmap[L->sn_struct[K][first_row+j]]-1),L->up_blocks[J][(bitmap[L->sn_struct[K][first_row+j]]-1)*(L->up_blocks_ld)[J]+(bitmap[L->sn_struct[K][ir+first_row]]-1)]);*/

	}
    /*
    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;
    */
    
    for (child = first_child[K]; child != -1; child = next_child[child]) {
      recursive_leftlooking_supernodal_update(J,child,
					      bitmap,dense_update_matrix,
					      A,L);
    }
  }

  /*
  else
    for(i=0;i<sn_up_size_father;i++)
      bitmap[L->sn_struct[J][i]]=0;
  */

}

static int
leftlooking_supernodal_front_factor(int sn,
				    int* bitmap,
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
  for(ip=0;ip<(L->sn_up_size)[sn];ip++) bitmap[(L->sn_struct)[sn][ip]] = ip;

  for(jp=0;jp<sn_size;jp++) {
    ind = &(A->rowind[A->colptr[ (L->sn_struct)[sn][jp] ]]);
    re  = &(A->taucs_values[A->colptr[ (L->sn_struct)[sn][jp] ]]); 
    for(ip=0;
	ip < A->colptr[ (L->sn_struct)[sn][jp] + 1 ] 
           - A->colptr[ (L->sn_struct)[sn][jp] ];
	ip++) {
      if (bitmap[ind[ip]] < sn_size)
	(L->sn_blocks)[sn][ (L->sn_blocks_ld)[sn]*jp + bitmap[ind[ip]]] =
	  taucs_add( (L->sn_blocks)[sn][ (L->sn_blocks_ld)[sn]*jp + bitmap[ind[ip]]] , re[ip] );
      else
	(L->up_blocks)[sn][ (L->up_blocks_ld)[sn]*jp + bitmap[ind[ip]] - sn_size] =
	taucs_add( (L->up_blocks)[sn][ (L->up_blocks_ld)[sn]*jp + bitmap[ind[ip]] - sn_size] , re[ip] );
    }
  }
  
  /* we use the BLAS through the Fortran interface */

  /* solving of lower triangular system for L */
  if (sn_size)
    taucs_potrf ("LOWER",
		 &sn_size,
		 (L->sn_blocks)[sn],&((L->sn_blocks_ld)[sn]),
		 &INFO);

  if (INFO) {
    taucs_printf("\t\tLL^T Factorization: Matrix is not positive definite.\n");
    taucs_printf("\t\t                    nonpositive pivot in column %d\n",
		 (L->sn_struct)[INFO-1]);
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
		(L->sn_blocks)[sn],&((L->sn_blocks_ld)[sn]),
		(L->up_blocks)[sn],&((L->up_blocks_ld)[sn]));

  return 0;
}

static int
recursive_leftlooking_supernodal_factor_llt(int sn,       /* this supernode */
					    int is_root,  /* is v the root? */
					    int* bitmap,
					    int* indmap,
					    taucs_ccs_matrix* A,
					    supernodal_factor_matrix* L)
{
  int  child;
  int  sn_size;
  int* first_child   = L->first_child;
  int* next_child    = L->next_child;
  taucs_datatype* dense_update_matrix = NULL;

  if (!is_root)
    sn_size = L->sn_size[sn];
  else
    sn_size = -1;

  if (!is_root) { 
    (L->sn_blocks   )[sn] = (L->up_blocks   )[sn] = NULL;
    if (L->sn_size[sn]) {
      (L->sn_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_size)[sn])*((L->sn_size)[sn]),
							    sizeof(taucs_datatype));
      if (!((L->sn_blocks)[sn])) return -1; /* the caller will free L */
    }
    (L->sn_blocks_ld)[sn] = (L->sn_size   )[sn];

    if (((L->sn_up_size)[sn] - (L->sn_size)[sn]) && (L->sn_size)[sn]) {
      (L->up_blocks   )[sn] = (taucs_datatype*)taucs_calloc(((L->sn_up_size)[sn]-(L->sn_size)[sn])
							    *((L->sn_size)[sn]),sizeof(taucs_datatype));
      if (!((L->up_blocks)[sn])) return -1; /* the caller will free L */
    }
    (L->up_blocks_ld)[sn] = (L->sn_up_size)[sn]-(L->sn_size)[sn];
  }

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    if (recursive_leftlooking_supernodal_factor_llt(child,
						    FALSE,
						    bitmap,
						    indmap,
						    A,L)
	== -1 ) {
      taucs_free(dense_update_matrix);
      return -1;
    }
    
    if (!is_root) {
      if (!dense_update_matrix) {
	dense_update_matrix = 
	  (taucs_datatype*) taucs_calloc((L->sn_up_size)[sn]*(L->sn_size)[sn],sizeof(taucs_datatype));
	if (!dense_update_matrix) return -1; /* caller will free L */
      }

      /* prepare the bitmap. Moved out of the recusive
	 update procedure 20/1/2003. Sivan and Elad */

      {
	int i;
	int J = sn;
	int sn_size_father = (L->sn_size)[J];
	int sn_up_size_father = (L->sn_up_size)[J];

	for(i=0;i<sn_size_father;i++)
	  bitmap[L->sn_struct[J][i]]=i+1;
	for(i=sn_size_father;i<sn_up_size_father;i++)
	  bitmap[L->sn_struct[J][i]] = i - sn_size_father + 1;
      }

      recursive_leftlooking_supernodal_update(sn,child,
					      bitmap,dense_update_matrix,
					      A,L);

      {
	int i;
	int J = sn;
	int sn_size_father = (L->sn_size)[J];
	int sn_up_size_father = (L->sn_up_size)[J];

	for(i=0;i<sn_size_father;i++)
	  bitmap[L->sn_struct[J][i]]=0;
	for(i=0;i<sn_up_size_father;i++)
	  bitmap[L->sn_struct[J][i]]=0;
      }

    }
  }
  taucs_free(dense_update_matrix);
  
  if(!is_root) {
    if (leftlooking_supernodal_front_factor(sn,
					    indmap,
					    A,
					    L)) {
      return -1; /* nonpositive pivot */
    }
  }

  return 0;
}


void* 
taucs_dtl(ccs_factor_llt_ll)(taucs_ccs_matrix* A)
{
  return taucs_dtl(ccs_factor_llt_ll_maxdepth)(A,0);
}

void* 
taucs_dtl(ccs_factor_llt_ll_maxdepth)(taucs_ccs_matrix* A,int max_depth)
{
  supernodal_factor_matrix* L;
  int* map;
  int *map2;
  double wtime, ctime;
  int fail;

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  L = multifrontal_supernodal_create();
  if (!L) return NULL;

  fail = taucs_ccs_symbolic_elimination(A,L,
					TRUE /* sort row indices */,
					max_depth);

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSymbolic Analysis            = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  map  = (int*)taucs_malloc((A->n+1)*sizeof(int));
  map2 = (int*)taucs_calloc((A->n+1),sizeof(int));

  if (fail == -1 || !map || !map2) {
    taucs_supernodal_factor_free(L);
    taucs_free(map2);
    taucs_free(map);
    return NULL;
  }

  wtime = taucs_wtime();
  ctime = taucs_ctime();

  if (recursive_leftlooking_supernodal_factor_llt((L->n_sn),  
						  TRUE, 
						  map2,
						  map,
						  A,L)
      == -1) {
    taucs_supernodal_factor_free(L);
    taucs_free(map);
    taucs_free(map2);
    return NULL;
  }

  wtime = taucs_wtime()-wtime;
  ctime = taucs_ctime()-ctime;
  taucs_printf("\t\tSupernodal Left-Looking LL^T = % 10.3f seconds (%.3f cpu)\n",
	       wtime,ctime);

  taucs_free(map);
  taucs_free(map2);

  return (void*) L;
}


/*************************************************************/
/* supernodal solve routines                                 */
/*************************************************************/

static void 
recursive_supernodal_solve_l(int sn,       /* this supernode */
			     int is_root,  /* is v the root? */
			     int* first_child, int* next_child,
			     int** sn_struct, int* sn_sizes, int* sn_up_sizes,
			     int* sn_blocks_ld,taucs_datatype* sn_blocks[],
			     int* up_blocks_ld,taucs_datatype* up_blocks[],
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
  int i;/*ip,j,jp omer*/

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    recursive_supernodal_solve_l(child,
				 FALSE,
				 first_child,next_child,
 				 sn_struct,sn_sizes,sn_up_sizes,
				 sn_blocks_ld,sn_blocks,
				 up_blocks_ld,up_blocks,
				 x,b,t);
  }

  if(!is_root) {

    sn_size = sn_sizes[sn];
    up_size = sn_up_sizes[sn] - sn_sizes[sn];

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
		  sn_blocks[sn],&(sn_blocks_ld[sn]),
		  xdense       ,&sn_size);
      
      if (up_size > 0 && sn_size > 0) {
	taucs_gemm ("No Conjugate","No Conjugate",
		    &up_size, &ione, &sn_size,
		    &taucs_one_const,
		    up_blocks[sn],&(up_blocks_ld[sn]),
		    xdense       ,&sn_size,
		    &taucs_zero_const,
		    bdense       ,&up_size);
      }
      
      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = xdense[i];
      for (i=0; i<up_size; i++)
	/*b[ sn_struct[ sn ][ sn_size + i ] ] -= bdense[i];*/
	b[ sn_struct[ sn ][ sn_size + i ] ] =
	  taucs_sub( b[ sn_struct[ sn ][ sn_size + i ] ] , bdense[i] );

#if 1
    }
#else
    } else if (sn_size > SOLVE_DENSE_CUTOFF) {

      xdense = t;
      bdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	xdense[i] = b[ sn_struct[ sn ][ i ] ];
      for (i=0; i<up_size; i++)
	bdense[i] = 0.0;
      
      for (jp=0; jp<sn_size; jp++) {
	xdense[jp] = xdense[jp] / sn_blocks[sn][ sn_blocks_ld[sn]*jp + jp];

	for (ip=jp+1; ip<sn_size; ip++) {
	  xdense[ip] -= xdense[jp] * sn_blocks[sn][ sn_blocks_ld[sn]*jp + ip];
	}
      }

      for (jp=0; jp<sn_size; jp++) {
	for (ip=0; ip<up_size; ip++) {
	  bdense[ip] += xdense[jp] * up_blocks[sn][ up_blocks_ld[sn]*jp + ip];
	}
      }

      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = xdense[i];
      for (i=0; i<up_size; i++)
	b[ sn_struct[ sn ][ sn_size + i ] ] -= bdense[i];
      
    } else {

      for (jp=0; jp<sn_size; jp++) {
	j = sn_struct[sn][jp];
	x[j] = b[j] / sn_blocks[sn][ sn_blocks_ld[sn]*jp + jp];
	for (ip=jp+1; ip<sn_size; ip++) {
	  i = sn_struct[sn][ip];
	  b[i] -= x[j] * sn_blocks[sn][ sn_blocks_ld[sn]*jp + ip];
	}

	for (ip=0; ip<up_size; ip++) {
	  i = sn_struct[sn][sn_size + ip];
	  b[i] -= x[j] * up_blocks[sn][ up_blocks_ld[sn]*jp + ip];
	}
      }

    }
#endif
  }
}

static void 
recursive_supernodal_solve_lt(int sn,       /* this supernode */
 			      int is_root,  /* is v the root? */
			      int* first_child, int* next_child,
			      int** sn_struct, int* sn_sizes, int* sn_up_sizes,
			      int* sn_blocks_ld,taucs_datatype* sn_blocks[],
			      int* up_blocks_ld,taucs_datatype* up_blocks[],
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
  int i;/*ip,j,jp omer*/

  if(!is_root) {

    sn_size = sn_sizes[sn];
    up_size = sn_up_sizes[sn]-sn_sizes[sn];
    
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
		     up_blocks[sn],&(up_blocks_ld[sn]),
		     xdense       ,&up_size,
		     &taucs_one_const,
		     bdense       ,&sn_size);
      
      taucs_trsm ("Left",
		  "Lower",
		  "Conjugate",
		  "No unit diagonal",
		  &sn_size,&ione,
		  &taucs_one_const,
		  sn_blocks[sn],&(sn_blocks_ld[sn]),
		  bdense       ,&sn_size);
      
      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = bdense[i];
#if 1
    }
#else    
    } else if (sn_size > SOLVE_DENSE_CUTOFF) {

      bdense = t;
      xdense = t + sn_size;
      
      for (i=0; i<sn_size; i++)
	bdense[i] = b[ sn_struct[ sn][ i ] ];
      for (i=0; i<up_size; i++)
	xdense[i] = x[ sn_struct[sn][sn_size+i] ];
      
      for (ip=sn_size-1; ip>=0; ip--) {
	for (jp=0; jp<up_size; jp++) {
	  bdense[ip] -= xdense[jp] * up_blocks[sn][ up_blocks_ld[sn]*ip + jp];
	}
      }

      for (ip=sn_size-1; ip>=0; ip--) {
	for (jp=sn_size-1; jp>ip; jp--) {
	  bdense[ip] -= bdense[jp] * sn_blocks[sn][ sn_blocks_ld[sn]*ip + jp];
	}
	bdense[ip] = bdense[ip] / sn_blocks[sn][ sn_blocks_ld[sn]*ip + ip];
      }

      for (i=0; i<sn_size; i++)
	x[ sn_struct[ sn][ i ] ]  = bdense[i];
    
    } else {

      for (ip=sn_size-1; ip>=0; ip--) {
	i = sn_struct[sn][ip];

	for (jp=0; jp<up_size; jp++) {
	  j = sn_struct[sn][sn_size + jp];
	  b[i] -= x[j] * up_blocks[sn][ up_blocks_ld[sn]*ip + jp];
	}

	for (jp=sn_size-1; jp>ip; jp--) {
	  j = sn_struct[sn][jp];
	  b[i] -= x[j] * sn_blocks[sn][ sn_blocks_ld[sn]*ip + jp];
	}
	x[i] = b[i] / sn_blocks[sn][ sn_blocks_ld[sn]*ip + ip];

      }

    }
#endif

  }

  for (child = first_child[sn]; child != -1; child = next_child[child]) {
    recursive_supernodal_solve_lt(child,
				  FALSE,
				  first_child,next_child,
				  sn_struct,sn_sizes,sn_up_sizes,
				  sn_blocks_ld,sn_blocks,
				  up_blocks_ld,up_blocks,
				  x,b,t);
  }
}


int 
taucs_dtl(supernodal_solve_llt)(void* vL, void* vx, void* vb)
{
  supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
  taucs_datatype* x = (taucs_datatype*) vx;
  taucs_datatype* b = (taucs_datatype*) vb;
  taucs_datatype* y;
  taucs_datatype* t; /* temporary vector */
  int     i;
  
  y = taucs_malloc((L->n) * sizeof(taucs_datatype));
  t = taucs_malloc((L->n) * sizeof(taucs_datatype));
  if (!y || !t) {
    taucs_free(y);
    taucs_free(t);
    taucs_printf("multifrontal_supernodal_solve_llt: out of memory\n");
    return -1;
  }

  for (i=0; i<L->n; i++) x[i] = b[i];

  recursive_supernodal_solve_l (L->n_sn,
				TRUE,  /* this is the root */
				L->first_child, L->next_child,
				L->sn_struct,L->sn_size,L->sn_up_size,
				L->sn_blocks_ld, L->sn_blocks,
				L->up_blocks_ld, L->up_blocks,
				y, x, t);

  recursive_supernodal_solve_lt(L->n_sn,
				TRUE,  /* this is the root */
				L->first_child, L->next_child,
				L->sn_struct,L->sn_size,L->sn_up_size,
				L->sn_blocks_ld, L->sn_blocks,
				L->up_blocks_ld, L->up_blocks,
				x, y, t);

  taucs_free(y);
  taucs_free(t);
    
  return 0;
}
#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*************************************************************/
/* generic interfaces to user-callable routines              */
/*************************************************************/

#ifdef TAUCS_CORE_GENERAL

cilk 
void* 
taucs_ccs_factor_llt_mf_maxdepth(taucs_ccs_matrix* A,int max_depth)
{
  void* p= NULL;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    p = spawn taucs_dccs_factor_llt_mf_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    p = spawn taucs_sccs_factor_llt_mf_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    p = spawn taucs_zccs_factor_llt_mf_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    p = spawn taucs_cccs_factor_llt_mf_maxdepth(A,max_depth);
#endif

  sync;
  return p;
}

void* 
taucs_ccs_factor_llt_ll_maxdepth(taucs_ccs_matrix* A,int max_depth)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt_ll_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt_ll_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt_ll_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt_ll_maxdepth(A,max_depth);
#endif

	assert(0);
  return NULL;
}

void* taucs_ccs_factor_llt_symbolic_maxdepth(taucs_ccs_matrix* A,int max_depth)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt_symbolic_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt_symbolic_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt_symbolic_maxdepth(A,max_depth);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt_symbolic_maxdepth(A,max_depth);
#endif
  
  assert(0);
  return NULL;
}

cilk
void* 
taucs_ccs_factor_llt_mf(taucs_ccs_matrix* A)
{
  void* p = NULL;
#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    p = spawn taucs_dccs_factor_llt_mf(A);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    p = spawn taucs_sccs_factor_llt_mf(A);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    p = spawn taucs_zccs_factor_llt_mf(A);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    p = spawn taucs_cccs_factor_llt_mf(A);
#endif

  sync;
  return p;
}

void* taucs_ccs_factor_llt_ll(taucs_ccs_matrix* A)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt_ll(A);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt_ll(A);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt_ll(A);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt_ll(A);
#endif
  
  assert(0);
  return NULL;
}

void* taucs_ccs_factor_llt_symbolic(taucs_ccs_matrix* A)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt_symbolic(A);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt_symbolic(A);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt_symbolic(A);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt_symbolic(A);
#endif
  
  assert(0);
  return NULL;
}

cilk
int taucs_ccs_factor_llt_numeric(taucs_ccs_matrix* A, void* L)
{
  int rc = TAUCS_ERROR;
#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    rc = spawn taucs_dccs_factor_llt_numeric(A,L);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    rc = spawn taucs_sccs_factor_llt_numeric(A,L);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    rc = spawn taucs_zccs_factor_llt_numeric(A,L);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    rc = spawn taucs_cccs_factor_llt_numeric(A,L);
#endif
  
  sync;
  return rc;
}


int taucs_supernodal_solve_llt(void* L, void* x, void* b)
{
	
#ifdef TAUCS_DOUBLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DOUBLE)
    return taucs_dsupernodal_solve_llt(L,x,b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SINGLE)
    return taucs_ssupernodal_solve_llt(L,x,b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DCOMPLEX)
    return taucs_zsupernodal_solve_llt(L,x,b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SCOMPLEX)
    return taucs_csupernodal_solve_llt(L,x,b);
#endif
  
  assert(0);
  return -1;
}

void taucs_supernodal_factor_free(void* L)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DOUBLE) {
    taucs_dsupernodal_factor_free(L);
    return;
  }
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SINGLE) {
    taucs_ssupernodal_factor_free(L);
    return;
  }
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DCOMPLEX) {
    taucs_zsupernodal_factor_free(L);
    return;
  }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SCOMPLEX) {
    taucs_csupernodal_factor_free(L);
    return;
  }
#endif
  
  assert(0);
}

void taucs_supernodal_factor_free_numeric(void* L)
{


#ifdef TAUCS_DOUBLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DOUBLE) {
    taucs_dsupernodal_factor_free_numeric(L);
    return;
  }
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SINGLE) {
    taucs_ssupernodal_factor_free_numeric(L);
    return;
  }
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DCOMPLEX) {
    taucs_zsupernodal_factor_free_numeric(L);
    return;
  }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SCOMPLEX) {
    taucs_csupernodal_factor_free_numeric(L);
    return;
  }
#endif
  
  assert(0);
}

taucs_ccs_matrix* 
taucs_supernodal_factor_to_ccs(void* L)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DOUBLE)
    return taucs_dsupernodal_factor_to_ccs(L);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SINGLE)
    return taucs_ssupernodal_factor_to_ccs(L);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DCOMPLEX)
    return taucs_zsupernodal_factor_to_ccs(L);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SCOMPLEX)
    return taucs_csupernodal_factor_to_ccs(L);
#endif
  
  assert(0);
  return NULL;
}

void* 
taucs_supernodal_factor_get_diag(void* L)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DOUBLE)
    return taucs_dsupernodal_factor_get_diag(L);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SINGLE)
    return taucs_ssupernodal_factor_get_diag(L);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_DCOMPLEX)
    return taucs_zsupernodal_factor_get_diag(L);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (((supernodal_factor_matrix*) L)->flags & TAUCS_SCOMPLEX)
    return taucs_csupernodal_factor_get_diag(L);
#endif
  
  assert(0);
  return NULL;
}

int
taucs_ccs_etree(taucs_ccs_matrix* A,
		int* parent,
		int* l_colcount,
		int* l_rowcount,
		int* l_nnz)
{
  int* prev_p;
  /*int* prev_nbr;omer*/
  int* level;
  /*int* first_descendant;omer*/
  int* l_cc;
  int* l_rc;
  int* wt;

  int  n = A->n;
  int  pprime;/*p,q,u omer*/
  int  ju;
  int* postorder;
  int* ipostorder;
  int  *first_child,*next_child;

  int i,j,k,ip,jp,kp;
  int nnz,jnnz;
  int* uf;
  int* rowptr;
  int* colind;
  int* rowcount;
  int* realroot;

  /* we need the row structures for the lower triangle */

  nnz = (A->colptr)[n];
  
  uf       = taucs_malloc(n     * sizeof(int));
  rowcount = taucs_malloc((n+1) * sizeof(int));
  rowptr   = taucs_malloc((n+1) * sizeof(int));
  colind   = taucs_malloc(nnz   * sizeof(int));

  if (!uf || !rowcount || !rowptr || !colind) {
    taucs_free(uf);
    taucs_free(rowcount);
    taucs_free(rowptr);
    taucs_free(colind);
    return -1;
  }

  for (i=0; i <=n; i++) rowcount[i] = 0;
  for (j=0; j < n; j++) {
    jnnz = (A->colptr)[j+1] - (A->colptr)[j];
    for (ip=0; ip<jnnz; ip++) {
      i = (A->rowind)[(A->colptr)[j] + ip];
      if (j < i) rowcount[i]++;
    }
  }

  ip = 0;
  for (i=0; i <= n; i++) {
    int next_ip = ip + rowcount[i];
    rowcount[i] = ip;
    rowptr  [i] = ip;
    ip = next_ip;
  }

  for (j=0; j < n; j++) {
    jnnz = (A->colptr)[j+1] - (A->colptr)[j];
    for (ip=0; ip<jnnz; ip++) {
      i = (A->rowind)[(A->colptr)[j] + ip];
      if (i==j) continue;
      assert( rowcount[i] < rowptr[i+1] );
      colind[ rowcount[i] ] = j;
      rowcount[i]++;
    }
  }

  /* now compute the etree */

  {
    int u,t,vroot;
    realroot = rowcount; /* reuse space */

    for (i=0; i<n; i++) {
      uf_makeset(uf,i);
      realroot[i] = i;
      parent[i] = n;
      vroot = i;
      for (kp=rowptr[i]; kp<rowptr[i+1]; kp++) {
	k=colind[kp];
	u = uf_find(uf,k);
	t = realroot[u];
	if (parent[t] == n && t != i) {
	  parent[t] = i;
	  vroot = uf_union(uf,vroot,u);
	  realroot[vroot] = i;
	}
      }
    }
  }

  taucs_free(colind);
  taucs_free(rowptr);
  taucs_free(rowcount);

  /* now only uf remains allocated */
  
  /* compute column counts */

  if (l_colcount || l_rowcount || l_nnz) {
    int* l_nz;
    int  tmp;
    int  u,p,q;

    first_child = taucs_malloc((n+1) * sizeof(int));
    next_child  = taucs_malloc((n+1) * sizeof(int));
    postorder   = taucs_malloc(n     * sizeof(int));
    ipostorder  = taucs_malloc(n     * sizeof(int));
    wt          = taucs_malloc(n     * sizeof(int));
    level       = taucs_malloc(n     * sizeof(int));
    prev_p      = taucs_malloc(n     * sizeof(int));

#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
    prev_nbr         = taucs_malloc(n     * sizeof(int));
    first_descendant = taucs_malloc(n     * sizeof(int));
#endif

    /* we allocate scratch vectors to avoid conditionals */
    /* in the inner loop.                                */

    if (l_colcount) l_cc = l_colcount;
    else            l_cc = (int*) taucs_malloc(n*sizeof(int));
    if (l_rowcount) l_rc = l_rowcount;
    else            l_rc = (int*) taucs_malloc(n*sizeof(int));
    if (l_nnz)      l_nz = l_nnz;
    else            l_nz = &tmp;


    if (!first_child || !next_child || !postorder
	|| !ipostorder || !wt || !level || !prev_p
	|| (!l_colcount && !l_cc) 
	|| (!l_rowcount && !l_rc) 
#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
	|| !prev_nbr || !first_descendant
#endif
	) {
      taucs_free(uf);

      if (!l_colcount) taucs_free(l_cc);
      if (!l_rowcount) taucs_free(l_rc);

      taucs_free(postorder);
      taucs_free(ipostorder);
      taucs_free(wt);
      taucs_free(level);
      taucs_free(prev_p);
      
#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
      taucs_free(prev_nbr);
      taucs_free(first_descendant);
#endif
      return -1;
    }

    /*for (j=0; j<n; j++) printf("parent[%d] = %d\n",j,parent[j]);*/

    /* compute the postorder */
    
    for (j=0; j<=n; j++) first_child[j] = -1;
    for (j=n-1; j>=0; j--) {
      next_child[j] = first_child[parent[j]];
      first_child[parent[j]] = j;
    }
    
    {
      int next = 0;
      recursive_postorder(n,first_child,next_child,
			  postorder,
			  ipostorder,&next);
    }
    
    /* sort by postorder of etree */
    /* compute level, fst_desc */

    tree_level(n,TRUE,first_child,next_child,
	       level,-1);
    
    for (u=0; u < n; u++) prev_p  [u] = -1;
    for (u=0; u < n; u++) l_rc    [u] =  1;
    for (u=0; u < n; u++) ordered_uf_makeset(uf,u);
    for (u=0; u < n; u++) {
      if (first_child[u] == -1)
	wt[u] = 1; /* leaves     */
      else
	wt[u] =  0; /* nonleaves */
    }

#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
    for (u=0; u < n; u++) prev_nbr[u] = -1;

    tree_first_descendant(n,TRUE,
			  first_child,next_child,ipostorder,
			  first_descendant);
#endif

    taucs_free(first_child);
    taucs_free(next_child);

    for (p=0; p<n; p++) {
      jp = postorder[p];
      if (parent[jp] != n) wt[parent[jp]] --;
      for (ip = (A->colptr)[jp]; ip < (A->colptr)[jp+1]; ip++) {
	ju = (A->rowind)[ip];
	u  = ipostorder[ju];
	if (ju==jp) continue; /* we only want proper neighbors */
#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
	if (first_descendant[jp] > prev_nbr[u]) {
#else
	if (1) {
#endif
	  wt[jp] ++;
	  pprime = prev_p[ju];
	  if (pprime == -1) 
	    l_rc[ju] += level[jp] - level[ju];
	  else {
	    q = ordered_uf_find(uf,pprime);
	    l_rc[ju] += level[jp] - level[q];
	    wt[q] --;
	  }
	  prev_p[ju] = jp;
	}
#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
	prev_nbr[u] = p;
#endif
      }
      if (parent[jp] != n) {
	if (!(ipostorder[parent[jp]] > ipostorder[jp])) {
	  printf("jp %d parent %d (ipo_j %d ipo_parent %d\n",
		 jp,parent[jp],ipostorder[jp],ipostorder[parent[jp]]);
	}
	assert(ipostorder[parent[jp]] > ipostorder[jp]);
	ordered_uf_union(uf,jp,parent[jp]);
      }
    }

    *l_nz = 0;
    for (u=0; u<n; u++) {
      l_cc[u] = wt[u];
      *l_nz += wt[u];
    }
    for (u=0; u<n; u++) {
      if (parent[u] != n) {
	l_cc[parent[u]] += l_cc[u];
	*l_nz += l_cc[u];
      }
    }

    /* free scrtach vectors                              */

    if (!l_colcount) taucs_free(l_cc);
    if (!l_rowcount) taucs_free(l_rc);

    /* free other data structures */

    taucs_free(postorder);
    taucs_free(ipostorder);
    taucs_free(wt);
    taucs_free(level);
    taucs_free(prev_p);
    
#ifdef GILBERT_NG_PEYTON_ANALYSIS_SUP
    taucs_free(prev_nbr);
    taucs_free(first_descendant);
#endif
  }

  taucs_free(uf);

  return 0;
}

int 
taucs_ccs_etree_liu(taucs_ccs_matrix* A,
		    int* parent,
		    int* l_colcount,
		    int* l_rowcount,
		    int* l_nnz)
{
  int n = A->n;
  int i,j,k,ip,kp;/*jp omer*/
  int nnz,jnnz;

  int* uf;
  int* rowptr;
  int* colind;

  int* rowcount;
  int* marker;
  int* realroot;

  int* l_cc;
  int* l_rc;

  /* we need the row structures for the lower triangle */

  nnz = (A->colptr)[n];
  
  uf       = taucs_malloc(n     * sizeof(int));
  rowcount = taucs_malloc((n+1) * sizeof(int));
  rowptr   = taucs_malloc((n+1) * sizeof(int));
  colind   = taucs_malloc(nnz   * sizeof(int));

  for (i=0; i <=n; i++) rowcount[i] = 0;

  for (j=0; j < n; j++) {
    
    jnnz = (A->colptr)[j+1] - (A->colptr)[j];

    for (ip=0; ip<jnnz; ip++) {
      i = (A->rowind)[(A->colptr)[j] + ip];
      if (j < i) rowcount[i]++;
    }
  }

  ip = 0;
  for (i=0; i <= n; i++) {
    int next_ip = ip + rowcount[i];
    rowcount[i] = ip;
    rowptr  [i] = ip;
    ip = next_ip;
  }

  for (j=0; j < n; j++) {
    jnnz = (A->colptr)[j+1] - (A->colptr)[j];

    for (ip=0; ip<jnnz; ip++) {
      i = (A->rowind)[(A->colptr)[j] + ip];
      if (i==j) continue;
      assert( rowcount[i] < rowptr[i+1] );
      colind[ rowcount[i] ] = j;
      rowcount[i]++;
    }
  }

  /* now compute the etree */

  {
    int u,t,vroot;
    realroot = rowcount; /* reuse space */

    for (i=0; i<n; i++) {
      uf_makeset(uf,i);
      realroot[i] = i;
      parent[i] = n;
      vroot = i;
      for (kp=rowptr[i]; kp<rowptr[i+1]; kp++) {
	k=colind[kp];
	u = uf_find(uf,k);
	t = realroot[u];
	if (parent[t] == n && t != i) {
	  parent[t] = i;
	  vroot = uf_union(uf,vroot,u);
	  realroot[vroot] = i;
	}
      }
    }
  }

  /* compute column counts */

  if (l_colcount || l_rowcount || l_nnz) {
    int* l_nz;
    int  tmp;

    /* we allocate scratch vectors to avoid conditionals */
    /* in the inner loop.                                */

    if (l_colcount) l_cc = l_colcount;
    else            l_cc = (int*) taucs_malloc(n*sizeof(int));
    if (l_rowcount) l_rc = l_rowcount;
    else            l_rc = (int*) taucs_malloc(n*sizeof(int));
    if (l_nnz)      l_nz = l_nnz;
    else            l_nz = &tmp;

    marker = rowcount; /* we reuse the space */
    
    for (j=0; j < n; j++) l_cc[j] = 1;
    *l_nz = n;
    
    for (i=0; i<n; i++) marker[i] = n; /* clear the array */
    
    for (i=0; i<n; i++) {
      l_rc[i] = 1;
      marker[ i ] = i;
      for (kp=rowptr[i]; kp<rowptr[i+1]; kp++) {
	k=colind[kp];
	j=k;
	while (marker[j] != i) {
	  l_cc[j]++;
	  l_rc[i]++;
	  (*l_nz)++;
	  marker[j] = i;
	  j = parent[j];
	}
      }
    }

    /* free scrtach vectors                              */

    if (!l_colcount) taucs_free(l_cc);
    if (!l_rowcount) taucs_free(l_rc);
  }

  taucs_free(colind);
  taucs_free(rowptr);
  taucs_free(rowcount);
  taucs_free(uf);

  return 0;
}


int
taucs_ccs_symbolic_elimination(taucs_ccs_matrix* A,
			       void* vL,
			       int do_order,
			       int max_depth
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

  int depth;

  L->n           = A->n;
  /* use calloc so we can deallocate unallocated entries */
  L->sn_struct   = (int**)taucs_calloc((A->n  ),sizeof(int*)); 
  L->sn_size     = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->sn_up_size  = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->first_child = (int*) taucs_malloc((A->n+1)*sizeof(int));
  L->next_child  = (int*) taucs_malloc((A->n+1)*sizeof(int));
  
  column_to_sn_map = (int*) taucs_malloc((A->n+1)*sizeof(int));
  map              = (int*) taucs_malloc((A->n+1)*sizeof(int));
  first_child      = (int*) taucs_malloc((A->n+1)*sizeof(int));
  next_child       = (int*) taucs_malloc((A->n+1)*sizeof(int));
  parent           = (int*) taucs_malloc((A->n+1)*sizeof(int));
  rowind           = (int*) taucs_malloc((A->n  )*sizeof(int));

  if (!(L->sn_struct) || !(L->sn_size) || !(L->sn_up_size) ||
      !(L->first_child) || !(L->next_child) || !column_to_sn_map
      || !map || !first_child || !next_child || !rowind || !parent) {
    taucs_free(parent);
    taucs_free(rowind);
    taucs_free(next_child);
    taucs_free(first_child);
    taucs_free(map);
    taucs_free(column_to_sn_map);
    taucs_free(L->next_child);
    taucs_free(L->first_child);
    taucs_free(L->sn_up_size);
    taucs_free(L->sn_size);
    taucs_free(L->sn_struct);
    L->sn_struct = NULL;
    L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
    return -1;
  }

  if (taucs_ccs_etree(A,parent,NULL,NULL,NULL) == -1) {
    taucs_free(parent);
    taucs_free(rowind);
    taucs_free(next_child);
    taucs_free(first_child);
    taucs_free(map);
    taucs_free(column_to_sn_map);
    taucs_free(L->next_child);
    taucs_free(L->first_child);
    taucs_free(L->sn_up_size);
    taucs_free(L->sn_size);
    taucs_free(L->sn_struct);
    L->sn_struct = NULL;
    L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
    return -1;
  }

  if (0) {
    double wtime;
    int *cc1,*cc2,*rc1,*rc2;
    int *p1;
    int nnz1,nnz2;

    cc1=(int*)taucs_malloc((A->n)*sizeof(int));
    cc2=(int*)taucs_malloc((A->n)*sizeof(int));
    rc1=(int*)taucs_malloc((A->n)*sizeof(int));
    rc2=(int*)taucs_malloc((A->n)*sizeof(int));
    p1 =(int*)taucs_malloc((A->n)*sizeof(int));

    wtime = taucs_wtime();
    taucs_ccs_etree_liu(A,parent,cc1,rc1,&nnz1);
    wtime = taucs_wtime() - wtime;
    printf("\t\t\tLiu Analysis = %.3f seconds\n",wtime);

    wtime = taucs_wtime();
    taucs_ccs_etree(A,p1,cc2,rc2,&nnz2);
    wtime = taucs_wtime() - wtime;
    printf("\t\t\tGNP Analysis = %.3f seconds\n",wtime);

    for (j=0; j<(A->n); j++) assert(parent[j]==p1[j]);
    for (j=0; j<(A->n); j++) {
      if (cc1[j]!=cc2[j]) printf("j=%d cc1=%d cc2=%d\n",j,cc1[j],cc2[j]);
      assert(cc1[j]==cc2[j]);
    }

    for (j=0; j<(A->n); j++) {
      if (rc1[j]!=rc2[j]) printf("j=%d rc1=%d rc2=%d\n",j,rc1[j],rc2[j]);
      assert(rc1[j]==rc2[j]);
    }

    if (nnz1!=nnz2) printf("nnz1=%d nnz2=%d\n",nnz1,nnz2);
    
    taucs_free(cc1); taucs_free(cc2); taucs_free(rc1); taucs_free(rc2);
  }

  for (j=0; j <= (A->n); j++) first_child[j] = -1;
  for (j = (A->n)-1; j >= 0; j--) {
    int p = parent[j];
    next_child[j] = first_child[p];
    first_child[p] = j;
  }

  /* let's compute the depth of the etree, to bail out if it is too deep */
  /* the whole thing will work better if we compute supernodal etrees    */

  {
    int next_depth_count;
    int this_depth_count;
    int child,i;

    int* this_depth = rowind; /* we alias rowind */
    int* next_depth = map;    /* and map         */
    int* tmp;

    this_depth[0] = A->n;
    this_depth_count = 1;
    next_depth_count = 0;
    depth = -1;

    while (this_depth_count) {
      for (i=0; i<this_depth_count; i++) {
	child = first_child[ this_depth[i] ];
	while (child != -1) {
	  next_depth[ next_depth_count ] = child;
	  next_depth_count++;
	  child = next_child[ child ];
	}
      }
      
      tmp = this_depth;
      this_depth = next_depth;
      next_depth = tmp;

      this_depth_count = next_depth_count;
      next_depth_count = 0;
      depth++;
    }
  }

  taucs_printf("\t\tElimination tree depth is %d\n",depth);

  if (max_depth && depth > max_depth) {
    taucs_printf("taucs_ccs_symbolic_elimination: etree depth %d, maximum allowed is %d\n",
		 depth, max_depth);
    taucs_free(parent);
    taucs_free(rowind);
    taucs_free(next_child);
    taucs_free(first_child);
    taucs_free(map);
    taucs_free(column_to_sn_map);
    taucs_free(L->next_child);
    taucs_free(L->first_child);
    taucs_free(L->sn_up_size);
    taucs_free(L->sn_size);
    taucs_free(L->sn_struct);
    L->sn_struct = NULL;
    L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
    return -1;
  }

  /*
  taucs_free(parent);
  ipostorder = (int*)taucs_malloc((A->n+1)*sizeof(int));
  */
  
  ipostorder = parent;
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

  L->n_sn = 0;
  for (j=0; j < (A->n); j++) map[j] = -1;
  for (j=0; j <= (A->n); j++) (L->first_child)[j] = (L->next_child)[j] = -1;
  
  if (recursive_symbolic_elimination(A->n,
				     A,
				     first_child,next_child,
				     &(L->n_sn),
				     L->sn_size,L->sn_up_size,L->sn_struct,
				     L->first_child,L->next_child,
				     rowind,
				     column_to_sn_map,
				     map,
				     do_order,ipostorder
				     )
      == -1) {
    for (j=0; j < (A->n); j++) taucs_free((L->sn_struct)[j]);

    taucs_free(parent);
    taucs_free(rowind);
    taucs_free(next_child);
    taucs_free(first_child);
    taucs_free(map);
    taucs_free(column_to_sn_map);
    taucs_free(L->next_child);
    taucs_free(L->first_child);
    taucs_free(L->sn_up_size);
    taucs_free(L->sn_size);
    taucs_free(L->sn_struct);
    L->sn_struct = NULL;
    L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
    return -1;
  }

  {
    double nnz   = 0.0;
    double flops = 0.0;
    int sn,i,colnnz;
    int bytes;

    bytes = 
      1*sizeof(char)                /* uplo             */
      + 2*sizeof(int)               /* n, n_sn          */
      + 3*(L->n_sn)*sizeof(int)     /* etree            */
      + 4*(L->n_sn)*sizeof(int)     /* block sizes, lda */
      + 1*(L->n_sn)*sizeof(int*)    /* row/col indices  */
      + 3*(L->n_sn)*sizeof(taucs_datatype*) /* actual blocks    */
      ;

    for (sn=0; sn<(L->n_sn); sn++) {
      bytes += (L->sn_up_size)[sn] * sizeof(int);    
      bytes += ((L->sn_size)[sn]*(L->sn_up_size)[sn]) * sizeof(taucs_datatype);

      for (i=0, colnnz = (L->sn_up_size)[sn]; 
	   i<(L->sn_size)[sn]; 
	   i++, colnnz--) {
	/* There was a bug here. I did not count muliply-adds in the
	   update part of the computation as 2 flops but one. */
	/*flops += ((double)(colnnz) - 1.0) * ((double)(colnnz) + 2.0) / 2.0;*/
	flops += 1.0 + ((double)(colnnz)) * ((double)(colnnz));
	nnz   += (double) (colnnz);
      }
    }
    taucs_printf("\t\tSymbolic Analysis of LL^T: %.2e nonzeros, %.2e flops, %.2e bytes in L\n",
		 nnz, flops, (float) bytes);
  }

  for (j=0; j < (A->n); j++) map[j] = -1;
  if (1)
  (void) recursive_amalgamate_supernodes((L->n_sn) - 1,
					 &(L->n_sn),
					 L->sn_size,L->sn_up_size,L->sn_struct,
					 L->first_child,L->next_child,
					 rowind,
					 column_to_sn_map,
					 map,
					 do_order,ipostorder
					 );


  {
    double nnz   = 0.0;
    double flops = 0.0;
    int sn,i,colnnz;
    int bytes;

    bytes = 
      1*sizeof(char)                /* uplo             */
      + 2*sizeof(int)               /* n, n_sn          */
      + 3*(L->n_sn)*sizeof(int)     /* etree            */
      + 4*(L->n_sn)*sizeof(int)     /* block sizes, lda */
      + 1*(L->n_sn)*sizeof(int*)    /* row/col indices  */
      + 3*(L->n_sn)*sizeof(taucs_datatype*) /* actual blocks    */
      ;

    for (sn=0; sn<(L->n_sn); sn++) {
      bytes += (L->sn_up_size)[sn] * sizeof(int);
      bytes += ((L->sn_size)[sn]*(L->sn_up_size)[sn]) * sizeof(taucs_datatype);

      for (i=0, colnnz = (L->sn_up_size)[sn]; 
	   i<(L->sn_size)[sn]; 
	   i++, colnnz--) {
	/* There was a bug here. I did not count muliply-adds in the
	   update part of the computation as 2 flops but one. */
	/*flops += ((double)(colnnz) - 1.0) * ((double)(colnnz) + 2.0) / 2.0;*/
	flops += 1.0 + ((double)(colnnz)) * ((double)(colnnz));
	nnz   += (double) (colnnz);
      }
    }
    taucs_printf("\t\tRelaxed  Analysis of LL^T: %.2e nonzeros, %.2e flops, %.2e bytes in L\n",
		 nnz, flops, (float) bytes);
  }

  /*
  {
    int i;
    printf("c2sn: ");
    for (i=0; i<A->n; i++) printf("%d ",column_to_sn_map[i]);
    printf("\n");
  }
  */
  
  taucs_free(parent);
  taucs_free(rowind);
  taucs_free(map);
  taucs_free(column_to_sn_map);
  taucs_free(next_child);
  taucs_free(first_child);

  L->sn_blocks_ld  = taucs_malloc((L->n_sn) * sizeof(int));
  L->sn_blocks     = taucs_calloc((L->n_sn), sizeof(taucs_datatype*)); /* so we can free before allocation */
  
  L->up_blocks_ld  = taucs_malloc((L->n_sn) * sizeof(int));
  L->up_blocks     = taucs_calloc((L->n_sn), sizeof(taucs_datatype*));

  if (!(L->sn_blocks_ld)
      || !(L->sn_blocks_ld)
      || !(L->sn_blocks)
      || !(L->up_blocks_ld)
      || !(L->up_blocks))
    return -1; /* the caller will free L */

  return 0;
}
#endif


/*************************************************************/
/* end of file                                               */
/*************************************************************/













