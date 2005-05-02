/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

typedef struct {
  int     length;
  int*    indices;
  int*    bitmap;

  taucs_datatype* values;
} spa;

#ifndef TAUCS_CORE_GENERAL

/*********************************************************/
/* sparse accumulator                                    */
/*********************************************************/

static spa* spa_create(int n)
{
  int i;
  spa* s;
  
  s = (spa*) taucs_malloc( sizeof(spa) );
  if ( !s ) return NULL;

  s->indices = (int*)    taucs_malloc( n * sizeof(int) );
  s->bitmap  = (int*)    taucs_malloc( n * sizeof(int) );
  s->values  = (taucs_datatype*) taucs_malloc( n * sizeof(taucs_datatype) );

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
  if (!s) return;

  taucs_free( s->indices );
  taucs_free( s->bitmap  );
  taucs_free( s->values  );
  taucs_free( s );
}


static void spa_set(spa* s, taucs_ccs_matrix* A, int j)
{
  int i, ip, next;
  taucs_datatype Aij;
  
  assert(j < A->n);

  next = 0;
  for (ip = (A->colptr)[j]; ip < (A->colptr)[j+1]; ip++) {
    i   = (A->rowind)[ip];
    Aij = (A->taucs_values)[ip];

    assert( i >= j ); /* A must be lower */
    
    (s->indices)[ next ] = i;
    (s->values) [ i    ] = Aij;
    (s->bitmap) [ i    ] = j;
    next++;
  }

  s->length = next;
}


static void spa_scale_add(spa* s, int j, taucs_ccs_matrix* A, int k, taucs_datatype alpha)
{
  int i, ip, next;
  taucs_datatype Aik;
  
  assert(k < A->n);

  /*
  printf("spa_scale_add: updating column %d with column %d\n",j,k);
  printf("spa_scale_add: colptr %d to %d-1\n",(A->colptr)[k],(A->colptr)[k+1]);
  */

  next = 0;
  for (ip = (A->colptr)[k]; ip < (A->colptr)[k+1]; ip++) {
    i   = (A->rowind)[ip];
    if (i < j) continue;
    Aik = (A->taucs_values)[ip];
    
    if ( (s->bitmap)[ i ] < j ) {
      /*printf("fill in (%d,%d)\n",i,j);*/
      (s->bitmap)[i] = j;
      (s->values)[i] = taucs_zero;
      (s->indices)[ s->length ] = i;
      (s->length)++;
    }

    /*(s->values)[ i ] += taucs_mul(alpha,Aik);*/

    (s->values)[ i ] = taucs_add((s->values)[ i ],
				 taucs_mul(alpha,Aik));

    /*printf("spa_scale_add: A(%d,%d) -= %lg * %lg ==> %lg\n",i,j,alpha,Aik,(s->values)[i]);*/
  }
}

/*********************************************************/
/* linked lists for rows                                 */
/*********************************************************/

static int*            rowlist;
static int*            rowlist_next;
static int*            rowlist_colind;
static taucs_datatype* rowlist_values;

static int     rowlist_freelist;
static int     rowlist_size;
static int     rowlist_next_expansion;

static int rowlist_create(int n)
{
  int i;

  rowlist_size           = 1000;
  rowlist_next_expansion = 1000;

  rowlist        = (int*) taucs_malloc( n * sizeof(int) );
  rowlist_next   = (int*) taucs_malloc( rowlist_size * sizeof(int) );
  rowlist_colind = (int*) taucs_malloc( rowlist_size * sizeof(int) );
  rowlist_values = (taucs_datatype*) taucs_malloc( rowlist_size * sizeof(taucs_datatype) );

  if (!rowlist || !rowlist_next | !rowlist_colind || !rowlist_values) {
    taucs_free(rowlist);
    taucs_free(rowlist_next);
    taucs_free(rowlist_colind);
    taucs_free(rowlist_values);
    rowlist = rowlist_next = rowlist_colind = NULL;
    rowlist_values = NULL;
    return -1;
  }

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

/* static void rowlist_freerow(int i){} */

static int rowlist_add(int i,int j,taucs_datatype v)
{
  int             l;
  int*            new_next;
  int*            new_colind;
  taucs_datatype* new_values;

  if (rowlist_freelist == -1) {
    int inc = rowlist_next_expansion;
    int ii;

    rowlist_next_expansion = (int) floor(1.25 * (double) rowlist_next_expansion);

    new_next   = (int*) taucs_realloc( rowlist_next,   (rowlist_size+inc) * sizeof(int) );
    if (!new_next) return -1;
    rowlist_next   = new_next;

    new_colind = (int*) taucs_realloc( rowlist_colind, (rowlist_size+inc) * sizeof(int) );
    if (!new_colind) return -1;
    rowlist_colind = new_colind;

    new_values = (taucs_datatype*) 
                            taucs_realloc(rowlist_values, 
				    (rowlist_size+inc) * sizeof(taucs_datatype) );
    if (!new_values) return -1;
    rowlist_values = new_values;

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


  return 0;
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

static taucs_datatype rowlist_getvalue(int l)
{
  return rowlist_values[ l ];
}

/*********************************************************/
/* Cholesky factorization                                */
/* This is a left-looking column-column code using       */
/* row lists. Can perform drop-tolerance incomplete      */
/* factorization with or without diagonal modification   */
/* to maintain rowsums.                                  */
/*********************************************************/

taucs_ccs_matrix* 
taucs_dtl(ccs_factor_llt)(taucs_ccs_matrix* A,double droptol, int modified)
{
  int            i,j,k,l,n,ip,next,Lnnz;
  taucs_datatype Lkj,pivot,v;
  double norm;
  spa*           s;
  taucs_ccs_matrix* L;
  taucs_datatype* dropped;
  int Aj_nnz;
  double flops = 0.0;

  if (!(A->flags & TAUCS_SYMMETRIC) && !(A->flags & TAUCS_HERMITIAN)) {

    taucs_printf("taucs_ccs_factor_llt: matrix must be symmetric\n");
    return NULL;
  }
  if (!(A->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_factor_llt: lower part must be represented\n");
    return NULL;
  }

  n = A->n;

  taucs_printf("taucs_ccs_factor_llt: starting n=%d droptol=%lf modified?=%d\n",
	     n,droptol,modified);

  L = taucs_dtl(ccs_create)(n,n,1000);
  if (!L) 
    return NULL;

  L->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

  Lnnz = 1000;
  next = 0;

  s = spa_create(n);
  i = rowlist_create(n);

  dropped = (taucs_datatype*) taucs_malloc( n * sizeof(taucs_datatype) );

  if (!s || i == -1 || !dropped) {
    taucs_ccs_free(L);
    spa_free(s);
    rowlist_free();
    taucs_free(dropped);
    return NULL;
  }

  for (i=0; i<n; i++) dropped[i] = taucs_zero;

  for (j=0; j<n; j++) {
    spa_set(s,A,j);

    for (l = rowlist_getfirst(j); 
	 l != -1; 
	 l = rowlist_getnext(l)) {
      k   = rowlist_getcolind(l);
      Lkj = rowlist_getvalue(l);
      /*spa_scale_add(s,j,L,k,taucs_neg(Lkj));*/ /* L_*j -= L_kj * L_*k */
      spa_scale_add(s,j,L,k,taucs_neg(taucs_conj(Lkj))); /* L_*j -= L_kj * L_*k */
    }

    /* we now add the j'th column of L to the taucs_ccs */
    
    if ( next+(s->length) > Lnnz ) {
      int*    rowind;
      taucs_datatype* values;
      int inc = max( (int) floor(1.25 * (double)Lnnz) , max( 8192, s->length ) );
      
      Lnnz += inc;

      rowind = (int*) taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      if (!rowind) {
	taucs_free(dropped);
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->rowind = rowind;

      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      if (!values) {
	taucs_free(dropped);
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->taucs_values = values;
    }

    (L->colptr)[j] = next;

    norm = 0.0;
    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      /* norm += v*v; */
      norm += taucs_re( taucs_mul(v,taucs_conj(v)) );
    }
    norm = sqrt(norm);

    Aj_nnz = (A->colptr)[j+1] - (A->colptr)[j];

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      if (i==j || taucs_abs(v) > droptol * norm || ip < Aj_nnz) {
	/* if (i==j || taucs_abs(v) > droptol * norm) { */
      } else {
	/*
	dropped[i] -= v;
	dropped[j] -= v;
	*/
	dropped[i] = taucs_sub( dropped[i], v );
	dropped[j] = taucs_sub( dropped[j], v );
      }
    }

    if (modified) 
      /*     pivot = taucs_sqrt( (s->values)[j] - dropped[j] );*/
      pivot = taucs_sqrt( taucs_sub ( (s->values)[j] , dropped[j] ) );
    else
      pivot = taucs_sqrt( (s->values)[j] );

#if 0
    taucs_printf("pivot=%.4e+%.4ei, sqrt=%.4e+%.4ei\n",
		 taucs_re( (s->values)[j] ),
		 taucs_im( (s->values)[j] ),
		 taucs_re( pivot ),
		 taucs_im( pivot ) );
#endif

    if (taucs_re(pivot) == 0.0 && taucs_im(pivot) == 0.0) {
      taucs_printf("taucs_ccs_factor_llt: zero pivot in column %d\n",j);
      taucs_printf("taucs_ccs_factor_llt: Ajj in spa = %lg dropped[j] = %lg Aj_nnz=%d\n",
		 (s->values)[j],dropped[j],Aj_nnz);
    } else if (taucs_abs(pivot) < 1e-12) {
      taucs_printf("taucs_ccs_factor_llt: small pivot in column %d (%le)\n",j,pivot);
    }

    /* we want Lii to be first in the compressed column */
    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];

      if (i==j) {
	/*if (modified) v = (s->values)[j] - dropped[j];*/
	if (modified) v = taucs_sub( (s->values)[j] , dropped[j] );

	/*v = v / pivot;*/
	v = taucs_div( v , pivot );

	(L->rowind)[next] = i;
	(L->taucs_values)[next] = v;
	next++;
	if (rowlist_add(i,j,v) == -1) {
	  taucs_free(dropped);
	  spa_free(s);
	  rowlist_free();
	  taucs_ccs_free(L);
	  return NULL;
	}
	break;
      }
    }

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      if (i==j) continue;

      /* if (modified && i == j) v = (s->values)[j] - dropped[j]; */

      if (i==j || taucs_abs(v) > droptol * norm || ip < Aj_nnz) {
	/* v = v / pivot; */
	v = taucs_div( v , pivot );

	(L->rowind)[next] = i;
	(L->taucs_values)[next] = v;
	next++;
	if (rowlist_add(i,j,v) == -1) {
	  taucs_free(dropped);
	  spa_free(s);
	  rowlist_free();
	  taucs_ccs_free(L);
	  return NULL;
	}
      }
    }

    (L->colptr)[j+1] = next;
    {
      double Lj_nnz = (double) ((L->colptr)[j+1] - (L->colptr)[j]);
      flops += 2.0 * Lj_nnz * Lj_nnz;
    }

    /* rowlist_free(j); */
  }

  (L->colptr)[n] = next;
  
  rowlist_free();
  spa_free(s);
  taucs_free(dropped);

  taucs_printf("taucs_ccs_factor_llt: done; nnz(L) = %d, flops=%.1le\n",(L->colptr)[n],flops);

  return L;
}

/***************** FACTOR LLT PARTIAL ********************/

/* 
 * Partial LL^T factorization. Factors the first p columns
 * and then updates, but does not factor, the trailing submatrix.
 * Designed for Shur-complement preconditioning.
 * 
 */

taucs_ccs_matrix* 
taucs_dtl(ccs_factor_llt_partial)(taucs_ccs_matrix* A, 
   			          int p)
{
  int            i,j,k,l,n,ip,next,Lnnz;
  taucs_datatype Lkj,pivot,v;
  spa*           s;
  taucs_ccs_matrix* L;
  int Aj_nnz;
  double flops = 0.0;

  if (!(A->flags & TAUCS_SYMMETRIC)) {
    taucs_printf("taucs_ccs_factor_llt_partial: matrix must be symmetric\n");
    return NULL;
  }
  if (!(A->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_factor_llt_partial: lower part must be represented\n");
    return NULL;
  }

  n = A->n;

  taucs_printf("taucs_ccs_factor_llt_partial: starting n=%d p=%d\n",n,p);

  L = taucs_dtl(ccs_create)(n,n,1000);
  if (!L) 
    return NULL;

  L->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

  Lnnz = 1000;
  next = 0;

  s = spa_create(n);
  i = rowlist_create(n);

  if (!s || i == -1) {
    taucs_ccs_free(L);
    spa_free(s);
    rowlist_free();
    return NULL;
  }

  for (j=0; j<p; j++) {
    spa_set(s,A,j);

    for (l = rowlist_getfirst(j); 
	 l != -1; 
	 l = rowlist_getnext(l)) {
      k   = rowlist_getcolind(l);
      Lkj = rowlist_getvalue(l);
      spa_scale_add(s,j,L,k,taucs_neg(Lkj)); /*  L_*j -= L_kj * L_*k  */
    }

    /* we now add the j'th column of L to the symccs */
    
    if ( next+(s->length) > Lnnz ) {
      int*    rowind;
      taucs_datatype* values;
      int inc = max( (int) floor(1.25 * (double)Lnnz) , max( 8192, s->length ) );
      /*int inc = max( 8192, s->length );*/
      
      Lnnz += inc;

      rowind = (int*) taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      if (!rowind) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->rowind = rowind;

      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      if (!values) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->taucs_values = values;

      /*
      rowind = (int*)    taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      assert( rowind && values );
      L->rowind = rowind;
      L->taucs_values = values;
      */
    }

    (L->colptr)[j] = next;

    Aj_nnz = (A->colptr)[j+1] - (A->colptr)[j]; 

    pivot = taucs_sqrt( (s->values)[j] );

    if (taucs_re(pivot) == 0.0 && taucs_im(pivot) == 0.0) {
      taucs_printf("taucs_ccs_factor_llt_partial: zero pivot in column %d\n",j);
    } else if (taucs_abs(pivot) < 1e-12) {
      taucs_printf("taucs_ccs_factor_llt_partial: small pivot in column %d (%le)\n",j,pivot);
    }

    /* we want Lii to be first in the compressed column */
    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];

      if (i==j) {
	/*v = v / pivot;*/
	v = taucs_div(v , pivot);

	(L->rowind)[next] = i;
	(L->taucs_values)[next] = v;
	next++;
	rowlist_add(i,j,v);
	break;
      }
    }

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      if (i==j) continue;

      /*v = v / pivot;*/
      v = taucs_div(v , pivot);

      (L->rowind)[next] = i;
      (L->taucs_values)[next] = v;
      next++;
      rowlist_add(i,j,v);
    }

    (L->colptr)[j+1] = next;

    {
      double Lj_nnz = (double) ((L->colptr)[j+1] - (L->colptr)[j]);
      flops += 2.0 * Lj_nnz * Lj_nnz;
    }
  }

  for (j=p; j<n; j++) {
    spa_set(s,A,j);

    /* we only apply updates from columns 0..p-1 */
    for (l = rowlist_getfirst(j); 
	 l != -1; 
	 l = rowlist_getnext(l)) {
      k   = rowlist_getcolind(l);
      Lkj = rowlist_getvalue(l);
      if (k >= p) continue; 
      spa_scale_add(s,j,L,k,taucs_neg(Lkj)); /*  L_*j -= L_kj * L_*k  */
    }

    /* we now add the j'th column of L to the symccs */
    
    if ( next+(s->length) > Lnnz ) {
      int*    rowind;
      taucs_datatype* values;
      int inc = max( (int) floor(1.25 * (double)Lnnz) , max( 8192, s->length ) );
      /*int inc = max( 8192, s->length );*/
      
      Lnnz += inc;

      rowind = (int*) taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      if (!rowind) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->rowind = rowind;

      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      if (!values) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->taucs_values = values;

      /*
      rowind = (int*)    taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      assert( rowind && values );
      L->rowind = rowind;
      L->taucs_values = values;
      */
    }

    (L->colptr)[j] = next;

    Aj_nnz = (A->colptr)[j+1] - (A->colptr)[j]; 

    /* we want Lii to be first in the compressed column */
    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];

      if (i==j) {
	(L->rowind)[next] = i;
	(L->taucs_values)[next] = v;
	next++;
	rowlist_add(i,j,v);
	break;
      }
    }

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      if (i==j) continue;
      (L->rowind)[next] = i;
      (L->taucs_values)[next] = v;
      next++;
      rowlist_add(i,j,v);
    }

    (L->colptr)[j+1] = next;

    /* not sure the flop count is correct. */
    {
      double Lj_nnz = (double) ((L->colptr)[j+1] - (L->colptr)[j]);
      flops += 2.0 * Lj_nnz * Lj_nnz;
    }
  }

  (L->colptr)[n] = next;
  
  rowlist_free();
  spa_free(s);

  taucs_printf("taucs_ccs_factor_llt_partial: done; nnz(L) = %d, flops=%.1le\n",(L->colptr)[n],flops);

  return L;
}

/*********************************************************/
/* LDL^T factorization                                   */
/* This is a left-looking column-column code using       */
/* row lists.                                            */
/*********************************************************/

taucs_ccs_matrix* 
taucs_dtl(ccs_factor_ldlt)(taucs_ccs_matrix* A)
{
  int            i,j,k,l,n,ip,next,Lnnz;
  taucs_datatype Lkj,pivot,v,Dkk;
  spa*           s;
  taucs_ccs_matrix* L;
  int Aj_nnz;
  double flops = 0.0;

  n = A->n;

  taucs_printf("taucs_ccs_factor_ldlt: starting n=%d\n",n);

  L = taucs_dtl(ccs_create)(n,n,1000);
  if (!L)
    return NULL;

  L->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

  Lnnz = 1000;
  next = 0;

  s = spa_create(n);
  i = rowlist_create(n);

  if (!s || i == -1) {
    taucs_ccs_free(L);
    spa_free(s);
    rowlist_free();
    return NULL;
  }

  for (j=0; j<n; j++) {
    spa_set(s,A,j);

    for (l = rowlist_getfirst(j); 
	 l != -1; 
	 l = rowlist_getnext(l)) {
      k   = rowlist_getcolind(l);
      Lkj = rowlist_getvalue(l);
      Dkk = (L->taucs_values)[ (L->colptr)[k] ];
      /*spa_scale_add(s,j,L,k,-Lkj*Dkk);*/ /* L_*j -= L_kj * L_*k */
      /*spa_scale_add(s,j,L,k,taucs_mul(taucs_neg(Lkj,Dkk));*/ /* L_*j -= L_kj * L_*k */
      spa_scale_add(s,j,L,k,
		    taucs_mul(taucs_neg(taucs_conj(Lkj)),Dkk)); /* L_*j -= L_kj * L_*k */
    }

    /* we now add the j'th column of L to the taucs_ccs */
    
    if ( next+(s->length) > Lnnz ) {
      int*    rowind;
      taucs_datatype* values;
      int inc = max( (int) floor(1.25 * (double)Lnnz) , max( 8192, s->length ) );
      /*int inc = max( 8192, s->length );*/
      
      Lnnz += inc;

      rowind = (int*) taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      if (!rowind) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->rowind = rowind;

      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      if (!values) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
      L->taucs_values = values;

      /*
      rowind = (int*)    taucs_realloc( L->rowind, Lnnz * sizeof(int) );
      values = (taucs_datatype*) taucs_realloc( L->taucs_values, Lnnz * sizeof(taucs_datatype) );
      assert( rowind && values );
      L->rowind = rowind;
      L->taucs_values = values;
      */
    }

    (L->colptr)[j] = next;

    Aj_nnz = (A->colptr)[j+1] - (A->colptr)[j]; 

    pivot = (s->values)[j]; 

    if (taucs_re(pivot) == 0.0 && taucs_im(pivot) == 0.0) {
      taucs_printf("ldlt: zero pivot in column %d\n",j);
      taucs_printf("ldlt: Ajj in spa = %lg Aj_nnz=%d\n",
		 (s->values)[j],Aj_nnz);
    }

#if 0
    if (taucs_abs(pivot) < 1e-4) {
      taucs_printf("taucs_ccs_factor_ldlt: small pivot in column %d: %.2e\n",
		   j,pivot);
    }
#endif

    /* we want Lii to be first in the compressed column */
    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];

      if (i==j) {
	/*
	printf(">>> %.8e + %.8ei / %.8e + %.8ei = ",
	       taucs_re(v),taucs_im(v),
	       taucs_re(pivot),taucs_im(pivot));
	*/

	/*v = v / pivot;*/
	v = taucs_div(v , pivot);

	/*printf("%.8e + %.8ei\n",taucs_re(v),taucs_im(v));*/

	(L->rowind)[next] = i;
	(L->taucs_values)[next] = pivot; /* we put D on the diagonal */
	next++;
	if (rowlist_add(i,j,v) == -1) {
	  spa_free(s);
	  rowlist_free();
	  taucs_ccs_free(L);
	  return NULL;
	}
	break;
      }
    }

    for (ip = 0; ip < s->length; ip++) {
      i = (s->indices)[ip];
      v = (s->values)[i];
      
      if (i==j) continue;

      /* v = v / pivot; */
      v = taucs_div(v , pivot);
      
      (L->rowind)[next] = i;
      (L->taucs_values)[next] = v;
      next++;
      if (rowlist_add(i,j,v) == -1) {
	spa_free(s);
	rowlist_free();
	taucs_ccs_free(L);
	return NULL;
      }
    }

    (L->colptr)[j+1] = next;
    {
      double Lj_nnz = (double) ((L->colptr)[j+1] - (L->colptr)[j]);
      flops += 2.0 * Lj_nnz * Lj_nnz;
    }

    /* rowlist_free(j); */
  }

  (L->colptr)[n] = next;
  
  rowlist_free();
  spa_free(s);

  taucs_printf("taucs_ccs_factor_ldlt: done; nnz(L) = %.2le, flops=%.2le\n",
	     (double) (L->colptr)[n],flops);

  return L;
}

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL
taucs_ccs_matrix* 
taucs_ccs_factor_llt(taucs_ccs_matrix* A,double droptol, int modified)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt(A,droptol,modified);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt(A,droptol,modified);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt(A,droptol,modified);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt(A,droptol,modified);
#endif
  
  assert(0);
  return NULL;
}

taucs_ccs_matrix* 
taucs_ccs_factor_llt_partial(taucs_ccs_matrix* A, int p)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_llt_partial(A,p);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_llt_partial(A,p);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_llt_partial(A,p);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_llt_partial(A,p);
#endif
  
  assert(0);
  return NULL;
}

taucs_ccs_matrix* 
taucs_ccs_factor_ldlt(taucs_ccs_matrix* A)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE)
    return taucs_dccs_factor_ldlt(A);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    return taucs_sccs_factor_ldlt(A);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_factor_ldlt(A);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_factor_ldlt(A);
#endif
  
  assert(0);
  return NULL;
}
#endif /*TAUCS_CORE_GENERAL*/

