/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/* File  : taucs_ccs_xxt.c                               */
/* Description: computes the Cholesky factor of A^-1     */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

#ifdef TAUCS_CORE_DOUBLE

typedef struct {
  int     length;
  int*    indices;
  int*    bitmap;
  double* values;
} spa;

/*********************************************************/
/* Returns the strictly upper part of A                  */
/*********************************************************/

static 
taucs_ccs_matrix*
ccs_syml_to_symu(taucs_ccs_matrix* A) {
  taucs_ccs_matrix* U;
  int n;
  int* temp;
  int i,j,ip;/*kp,k,jp omer*/
  double v;

  n = A->n;

  temp = (int*) taucs_malloc(n * sizeof(int));
  if (!temp) return NULL;

  U = taucs_dtl(ccs_create)(n, n, (A->colptr)[n] - n);
  if (!U) {
    taucs_free(temp);
    return NULL;
  }

  U->flags = TAUCS_SYMMETRIC | TAUCS_UPPER;

  for (j=0; j<=n; j++) (U->colptr)[j] = 0;
  for (j=0; j<n; j++)  temp[j] = 0;

  for (j=0; j<n; j++) {
    for (ip=(A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (i!=j) temp[i]++;
    }
  }

  for (j=1; j<=n; j++) (U->colptr)[j] = (U->colptr)[j-1] + temp[j-1];
  for (j=0; j< n; j++) temp[j] = (U->colptr)[j];
  
  for (j=0; j<n; j++) {
    for (ip=(A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      v = (A->taucs_values)[ip];
      if (i!=j) {
	(U->rowind)[ temp[i] ] = j;
	(U->taucs_values)[ temp[i] ] = v;
	temp[i]++;
      }
    }
  }

  /*
  taucs_ccs_write_ijv(A,"AA.ijv");
  taucs_ccs_write_ijv(U,"UU.ijv");
  */

  assert((U->colptr)[n] == (A->colptr)[n] - n);

  return U;
}

/*********************************************************/
/*                                                       */
/*********************************************************/

static spa* spa_create(int n)
{
  int i;
  spa* s;
  
  s = (spa*) taucs_malloc( sizeof(spa) );
  if ( !s ) return NULL;

  s->indices = (int*)    taucs_malloc( n * sizeof(int) );
  s->bitmap  = (int*)    taucs_malloc( n * sizeof(int) );
  s->values  = (double*) taucs_malloc( n * sizeof(double) );

  if ( !(s->indices) || !(s->values) || !(s->bitmap) ) {
    taucs_printf("chol: cannot create spa\n");
    taucs_free( s->indices );
    taucs_free( s->bitmap  );
    taucs_free( s->values  );
    taucs_free( s );
    return NULL;
  }

  s->length = 0;

  for (i=0; i<n; i++) (s->bitmap)[i] = -1;
  
  return s;
}

static void spa_free(spa* s)
{
  taucs_free( s->indices );
  taucs_free( s->values  );
  taucs_free( s );
}

static void spa_set_lu(spa* s, taucs_ccs_matrix* L, taucs_ccs_matrix* U, int j)
{
  int i, ip, next;
  double Aij;
  
  assert(j < L->n);

  next = 0;
  for (ip = (U->colptr)[j]; ip < (U->colptr)[j+1]; ip++) {
    i   = (U->rowind)[ip];
    Aij = (U->taucs_values)[ip];

    assert( i < j ); /* U must be strictly upper */
    
    (s->indices)[ next ] = i;
    (s->values) [ i    ] = Aij;
    (s->bitmap) [ i    ] = j;
    next++;
  }
  for (ip = (L->colptr)[j]; ip < (L->colptr)[j+1]; ip++) {
    i   = (L->rowind)[ip];
    Aij = (L->taucs_values)[ip];

    assert( i >= j ); /* A must be lower */
    
    (s->indices)[ next ] = i;
    (s->values) [ i    ] = Aij;
    (s->bitmap) [ i    ] = j;
    next++;
  }

  s->length = next;
}

static void spa_scale_add(spa* s, int j, taucs_ccs_matrix* A, int k, double alpha)
{
  int i, ip, next;
  double Aik;
  
  assert(k < A->n);

  /*
  printf("spa_scale_add: updating column %d with column %d\n",j,k);
  printf("spa_scale_add: colptr %d to %d-1\n",(A->colptr)[k],(A->colptr)[k+1]);
  */

  next = 0;
  for (ip = (A->colptr)[k]; ip < (A->colptr)[k+1]; ip++) {
    i   = (A->rowind)[ip];
    /*if (i < j) continue;*/
    Aik = (A->taucs_values)[ip];

    if ( (s->bitmap)[ i ] < j ) {
      /*printf("fill in (%d,%d)\n",i,j);*/
      (s->bitmap)[i] = j;
      (s->values)[i] = 0.0;
      (s->indices)[ s->length ] = i;
      (s->length)++;
    }

    (s->values)[ i ] += alpha * Aik;

    /*printf("spa_scale_add: A(%d,%d) -= %lg * %lg ==> %lg\n",i,j,alpha,Aik,(s->values)[i]);*/
  }
}
		    
static double spa_dot(spa* s, int j, taucs_ccs_matrix* A, int k)
{
  int i, ip;
  double Aik;
  double x = 0.0;
  
  assert(k < A->n);

  /*
  printf("spa_dot: updating column %d with column %d\n",j,k);
  printf("spa_dot: colptr %d to %d-1\n",(A->colptr)[k],(A->colptr)[k+1]);
  */

  for (ip = (A->colptr)[k]; ip < (A->colptr)[k+1]; ip++) {
    i   = (A->rowind)[ip];
    Aik = (A->taucs_values)[ip];
    

    if ( (s->bitmap)[ i ] == j ) {
      /*printf("j=%d, i=%d k=%d ::: %lg, %lg\n",j,i,k,Aik,(s->values)[ i ]);*/
      x += Aik * (s->values)[ i ];
    } else {
      /*printf("@@@ j=%d, i=%d k=%d ::: %lg, %lg\n",j,i,k,Aik,(s->values)[ i ]);*/
    }
  }

  return x;
}
		    
static double spa_A_norm(spa* s, int j, taucs_ccs_matrix* A)
{
  int i, ip, k, kp;
  double Aik;
  double x = 0.0;
  
  assert(A->flags | TAUCS_SYMMETRIC);
  assert(A->flags | TAUCS_LOWER);

  /*
  printf("spa_scale_add: updating column %d with column %d\n",j,k);
  printf("spa_scale_add: colptr %d to %d-1\n",(A->colptr)[k],(A->colptr)[k+1]);
  */

  for (kp=0; kp<s->length; kp++) {
    k = (s->indices)[kp];
    
    for (ip = (A->colptr)[k]; ip < (A->colptr)[k+1]; ip++) {
      i   = (A->rowind)[ip];
      Aik = (A->taucs_values)[ip];
      
      if ( (s->bitmap)[ i ] == j ) {
	/*printf("j=%d, i=%d k=%d ::: %lg, %lg\n",j,i,k,Aik,(s->values)[ i ]);*/
	if (i == k)
	  x += (s->values)[k] * Aik * (s->values)[ i ];
	else 
	  x += 2.0 * (s->values)[k] * Aik * (s->values)[ i ];
      }
    }
  }

  return x;
}
		    
/*********************************************************/
/*                                                       */
/*********************************************************/

int*    rowlist;
int*    rowlist_next;
int*    rowlist_colind;
double* rowlist_values;

int     rowlist_freelist;
int     rowlist_size;

static int rowlist_create(int n)
{
  int i;

  rowlist = (int*) taucs_malloc( n * sizeof(int) );

  rowlist_size    = 1000;
  rowlist_next    = (int*)    taucs_malloc( rowlist_size * sizeof(int) );
  rowlist_colind = (int*)    taucs_malloc( rowlist_size * sizeof(int) );
  rowlist_values  = (double*) taucs_malloc( rowlist_size * sizeof(double) );

  for (i=0; i<n; i++) rowlist[i] = -1; /* no list yet for row i */

  /* free list */
  rowlist_freelist = 0; 
  for (i=0; i<rowlist_size-1; i++) rowlist_next[i] = i+1; 
  rowlist_next[rowlist_size-1] = -1;
				   
  return 0;
}

static void rowlist_free()
{
  taucs_free(rowlist);
  taucs_free(rowlist_next);
  taucs_free(rowlist_colind);
  taucs_free(rowlist_values);
}

/*static void rowlist_freerow(int i){}*/

static void rowlist_add(int i,int j,double v)
{
  int l;

  if (rowlist_freelist == -1) {
    int inc = 1000;
    int ii;

    rowlist_next   = (int*)    taucs_realloc( rowlist_next,   (rowlist_size+inc) * sizeof(int) );
    rowlist_colind = (int*)    taucs_realloc( rowlist_colind, (rowlist_size+inc) * sizeof(int) );
    rowlist_values = (double*) taucs_realloc( rowlist_values, (rowlist_size+inc) * sizeof(double) );

    rowlist_freelist = rowlist_size;
    for (ii=rowlist_size; ii<rowlist_size+inc-1; ii++)
      rowlist_next[ii] = ii+1;
    rowlist_next[ rowlist_size+inc-1 ] = -1;

    rowlist_size    += inc;
  }

  l = rowlist_freelist;
  rowlist_freelist = rowlist_next[ rowlist_freelist ];

  rowlist_next  [ l ] = rowlist[ i ];
  rowlist_colind[ l ] = j;
  rowlist_values[ l ] = v;
  
  rowlist[ i ] = l;
}

static int rowlist_getfirst(int i)
{
  return rowlist[ i ];
}

static int rowlist_getnext(int l)
{
  return rowlist_next[ l ];
}

static int rowlist_getcolind(int l)
{
  return rowlist_colind[ l ];
}

/*
static double rowlist_getvalue(int l)
{
  return rowlist_values[ l ];
}
*/

/*********************************************************/
/* Inverse Cholesky factorization                        */
/*********************************************************/

taucs_ccs_matrix* 
taucs_ccs_factor_xxt(taucs_ccs_matrix* A)
{
  int            i,j,k,l,n,ip,next,Lnnz;
  double v;/*Lkj,pivot,norm omer*/
  spa*           s;
  spa*           Aej;
  taucs_ccs_matrix* L;
  taucs_ccs_matrix* U;
  /*int Aj_nnz;omer*/
  /*double flops = 0.0;*/
  double x;
  int* bitmap;

  if (!(A->flags & TAUCS_SYMMETRIC)) {
    taucs_printf("taucs_ccs_factor_xxt: matrix must be symmetric\n");
    return NULL;
  }
  if (!(A->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_factor_xxt: lower part must be represented\n");
    return NULL;
  }

  if (!(A->flags & TAUCS_DOUBLE)) {
    taucs_printf("taucs_ccs_factor_xxt: only works for double-precision real matrices\n");
    return NULL;
  }

  n = A->n;

  taucs_printf("taucs_ccs_factor_xxt: starting n=%d\n",n);

  bitmap = (int*) taucs_malloc(n * sizeof(int));
  if (!bitmap) return NULL;
  for (i=0; i<n; i++) bitmap[i] = -1;

  U = ccs_syml_to_symu(A);


  L = taucs_dtl(ccs_create)(n,n,1000);
  /*  L->flags = TAUCS_TRIANGULAR | TAUCS_LOWER; */
  L->flags = 0;

  Lnnz = 1000;
  next = 0;

  s   = spa_create(n);
  Aej = spa_create(n);
  rowlist_create(n);

  for (j=0; j<n; j++) {

    /* set the spa to ej */

    s->length = 1;
    (s->values)[j] = 1.0;
    (s->bitmap)[j] = j;
    (s->indices)[0] = j;

    /* compute A*ej, get both upper and lower parts! */

    spa_set_lu(Aej,A,U,j);

    /*for (k=0; k<j; k++) {*/

    for (ip=0; ip<Aej->length; ip++) {
      i = (Aej->indices)[ip];
      
      for (l = rowlist_getfirst(i); 
	   l != -1; 
	   l = rowlist_getnext(l)) {
	k   = rowlist_getcolind(l);
	
	if (bitmap[k] == j) continue;
	bitmap[k] = j;

	/* inner product of column k of X with A*ej */
      
	x = spa_dot(Aej,j,L,k);
	if (x != 0.0) {
	  /*printf("adding column %d to e_%d, before=%d\n",k,j,s->length);*/
	  spa_scale_add(s,j,L,k,-x); /* L_*j -= x * L_*k  */
	}
      }
    }

    /* normalize the column to unit A-norm */

    x = sqrt(spa_A_norm(s,j,A));
    /*printf("A-norm of column %d = %lg\n",j,x);*/

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      (s->values)[i] /= x;
    }

    /* we now add the j'th column of L to the taucs_ccs */
    
    if ( next+(s->length) > Lnnz ) {
      int*    rowind;
      double* values;
      int inc = max( 8192, s->length );
      
      Lnnz += inc;

      rowind = (int*)    taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      values = (double*) taucs_realloc( L->taucs_values, Lnnz * sizeof(double) );
      /* check for errors */
      assert( rowind && values );
      L->rowind = rowind;
      L->taucs_values = values;
    }

    (L->colptr)[j] = next;

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      (L->rowind)[next] = i;
      (L->taucs_values)[next] = v;
      next++;
      rowlist_add(i,j,v);
    }

    (L->colptr)[j+1] = next;
  }

  (L->colptr)[n] = next;
  
  taucs_free(bitmap);
  rowlist_free();
  spa_free(Aej);
  spa_free(s);
  taucs_ccs_free(U);

  taucs_printf("taucs_ccs_factor_xxt: done; nnz(L) = %d\n",(L->colptr)[n]);

  return L;
}

/*********************************************************/
/* XXT Solve                                             */
/*********************************************************/

int
taucs_ccs_solve_xxt(void* vX, double* x, double* b)
{
  taucs_ccs_matrix* X = (taucs_ccs_matrix*) vX;
  int n;
  int i,j,ip;
  double v;
  double* y;

  if (!(X->flags & TAUCS_TRIANGULAR)
      || !(X->flags & TAUCS_LOWER)
      || !(X->flags & TAUCS_DOUBLE)
      ) {
    taucs_printf("taucs_ccs_solve_xxt: matrix must be lower triangular double-precision real\n");
    return 0;
  }

  n = X->n;

  y = (double*) taucs_malloc(n * sizeof(double));
  if (!y) return -1;

  /* multiply by X' */

  for (j=0; j<n; j++) {
    y[j] = 0.0;

    for (ip=(X->colptr)[j]; ip<(X->colptr)[j+1]; ip++) {
      i = (X->rowind)[ip];
      v = (X->taucs_values)[ip];
      y[j] += v*b[i];
    }
  }

  for (i=0; i<n; i++) x[i] = 0.0;

  /* multiply by X */

  for (j=0; j<n; j++) {
    for (ip=(X->colptr)[j]; ip<(X->colptr)[j+1]; ip++) {
      i = (X->rowind)[ip];
      v = (X->taucs_values)[ip];
      x[i] += v*y[j];
    }
  }

  taucs_free(y);

  return 0;
}

#endif /* TAUCS_CORE_DOUBLE */
