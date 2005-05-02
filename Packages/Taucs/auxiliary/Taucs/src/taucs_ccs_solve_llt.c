/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "taucs.h"


/*********************************************************/
/*                                                       */
/*********************************************************/


#ifndef TAUCS_CORE_GENERAL
int 
taucs_dtl(ccs_solve_llt)(void* vL, taucs_datatype* x, taucs_datatype* b)
{
  taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

  int n;
  int i,j;
  int ip,jp;
  taucs_datatype  Aij, Ajj, Aii;
  taucs_datatype* y;

  if (!(L->flags & TAUCS_TRIANGULAR)) {
    taucs_printf("taucs_ccs_solve_llt: factor matrix must be triangular\n");
    return -1;
  }
  if (!(L->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_solve_llt: lower part must be represented\n");
    return -1;
  }

  n = L->n;

  y = (taucs_datatype*) taucs_malloc( n * sizeof(taucs_datatype) );
  if (!y) return -1;

  for (i=0; i<n; i++) x[i] = b[i];

  /* Solve L y = b = x  */

  for (j=0; j<n; j++) {

    /* we put diagonal elements first */
    ip = (L->colptr)[j];
    i = (L->rowind)[ip];
    assert (i==j);
    Ajj = (L->taucs_values)[ip];
    
    /* 
    for (ip = (L->colptr)[j]; ip < (L->colptr)[j+1]; ip++) {
      i = (L->rowind)[ip];
      if (i==j) {
	Ajj = (L->taucs_values)[ip];
	break;
      }
    }
    */

    /*y[j] = x[j] / Ajj;*/
    y[j] = taucs_div( x[j] , Ajj );

    for (ip = (L->colptr)[j] + 1; ip < (L->colptr)[j+1]; ip++) {
      i = (L->rowind)[ip];
      Aij = (L->taucs_values)[ip];
      /*x[i] -= y[j]*Aij;*/
      x[i] = taucs_sub( x[i], taucs_mul( y[j],Aij ));
    }

    /*
    for (ip = (L->colptr)[j]; ip < (L->colptr)[j+1]; ip++) {
      i = (L->rowind)[ip];
      if (i != j) {
	Aij = (L->taucs_values)[ip];
	x[i] -= y[j]*Aij;
      }
    }
    */
  }
  
  /* Solve L^T x = y */

  for (i=n-1; i>=0; i--) {

    for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++) {
      j = (L->rowind)[jp];
      Aij = taucs_conj( (L->taucs_values)[jp] );
      /*y[i] -= x[j]*Aij;*/
      y[i] = taucs_sub( y[i], taucs_mul( x[j],Aij ));
    }
    /*
    for (jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++) {
      j = (L->rowind)[jp];
      if (i != j) {
	Aij = (L->taucs_values)[jp];
	y[i] -= x[j]*Aij;
      }
    }
    */

    jp = (L->colptr)[i];
    j = (L->rowind)[jp];
    Aii = (L->taucs_values)[jp]; 

    /*
    for (jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++) {
      j = (L->rowind)[jp];
      if (i==j) {
	Aii = (L->taucs_values)[jp];
	break;
      }
    }
    */

    /*x[i] = y[i] / Aii;*/
    x[i] = taucs_div( y[i] , Aii );

  }

  taucs_free(y);

  return 0;
}

/***************** SOLVE LLT PARTIAL ********************/

int 
taucs_dtl(ccs_solve_schur)(taucs_ccs_matrix* L,
			   taucs_ccs_matrix* schur_comp,
			   int    (*schur_precond_fn)(void*,void* x,void* b),
			   void*  schur_precond_args,
			   int    maxits,
			   double convratio,
			   taucs_datatype* x, taucs_datatype* b)
{
  int n;
  int i,j;
  int ip,jp;
  taucs_datatype  Aij, Ajj, Aii;
  taucs_datatype* y;

  int p;

  if (!(L->flags & TAUCS_TRIANGULAR)) {
    taucs_printf("taucs_ccs_solve_llt: factor matrix must be triangular\n");
    return -1;
  }
  if (!(L->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_solve_llt: lower part must be represented\n");
    return -1;
  }

  n = L->n;
  p = n - (schur_comp->n);

  y = (taucs_datatype*) taucs_malloc( n * sizeof(taucs_datatype) );
  if (!y) return -1;

  for (i=0; i<n; i++) x[i] = b[i];

  /* Solve L y = b = x  */

  for (j=0; j<p; j++) {

    /* we put diagonal elements first */
    ip = (L->colptr)[j];
    i = (L->rowind)[ip];
    assert (i==j);
    Ajj = (L->taucs_values)[ip];
    
    /*y[j] = x[j] / Ajj;*/
    y[j] = taucs_div( x[j] , Ajj );

    for (ip = (L->colptr)[j] + 1; ip < (L->colptr)[j+1]; ip++) {
      i = (L->rowind)[ip];
      Aij = (L->taucs_values)[ip];
      /*x[i] -= y[j]*Aij;*/
      x[i] = taucs_sub( x[i], taucs_mul( y[j],Aij ));
    }
  }

  /* 
     now y_1 is computed, L_11 y_1 = b_1, 
     x_2 holds (b_2 - L_21 y_1).
     move y_2 <- x_2.
  */

  for (i=p; i<n; i++) y[i] = x[i];

  /* Now solve x_2 <- (A_22 - L_21 L_21^T)^-1 y_2 */ 
  
  /*taucs_printf("symccs_solve_schur: calling CG on Schur complement\n");*/
  /* sivan: removed for testing the complex codes */
  assert(0);
#if 0
  taucs_conjugate_gradients (schur_comp,
			     schur_precond_fn,
			     schur_precond_args,
			     x+p,          /* this is x_2 */
			     y+p,          /* this is y_2 */
			     maxits,       /* itermax */
			     convratio     /* conv tolerance */
			     );
#endif
  /*taucs_printf("taucs_ccs_solve_llt_partial: CG on Schur complement returned\n");*/
  
  /* Now we have x_2, solve L_11^T x_1 = y_1 - L_21^T x_2 */

  for (i=p-1; i>=0; i--) {

    for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++) {
      j = (L->rowind)[jp];
      Aij = (L->taucs_values)[jp];
      /*y[i] -= x[j]*Aij;*/
      y[i] = taucs_sub( y[i], taucs_mul( x[j],Aij ));
    }

    jp = (L->colptr)[i];
    j = (L->rowind)[jp];
    Aii = (L->taucs_values)[jp];

    /*x[i] = y[i] / Aii;*/
    x[i] = taucs_div( y[i] , Aii );

  }

  taucs_free(y);

  return 0;
}

/*********************************************************/
/* LDL^T solve                                           */
/*********************************************************/

int 
taucs_dtl(ccs_solve_ldlt)(void* vL, taucs_datatype* x, taucs_datatype* b)
{
  taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

  int n;
  int i,j;
  int ip,jp;
  taucs_datatype  Ajj = taucs_zero_const; /* just to suppress the warning */
  taucs_datatype  Aij = taucs_zero_const; /* just to suppress the warning */
  taucs_datatype* y;

  /* taucs_printf("taucs_ccs_solve_ldlt: starting\n"); */

  if (!(L->flags & TAUCS_TRIANGULAR)) {
    taucs_printf("taucs_ccs_solve_ldlt: factor matrix must be triangular\n");
    return -1;
  }
  if (!(L->flags & TAUCS_LOWER)) {
    taucs_printf("taucs_ccs_solve_ldlt: lower part must be represented\n");
    return -1;
  }

  n = L->n;

  y = (taucs_datatype*) taucs_malloc( n * sizeof(taucs_datatype) );
  if (!y) return -1;

  for (i=0; i<n; i++) x[i] = b[i];

  /* Solve L y = b = x  */

  /* taucs_printf("taucs_ccs_solve_ldlt: solving L y = b\n"); */

  for (j=0; j<n; j++) {

#if 0
    /* we put diagonal elements first */
    ip = (L->colptr)[j];
    i = (L->rowind)[ip];
    assert (i==j);
    /*Ajj = 1.0;*/
    Ajj = taucs_one;
    
    /*y[j] = x[j] / Ajj;*/
    y[j] = taucs_div( x[j] , Ajj );
#else
    y[j] = x[j];
#endif

    if (taucs_isnan(y[j]) || taucs_isinf(y[j])) {
      taucs_printf("taucs_ccs_solve_ldlt: inf/nan in column %d (L); %e+%ei / %e+%ei\n",
		   j,
		   taucs_re(x[j]),taucs_im(x[j]),
		   taucs_re(Ajj ),taucs_im(Ajj ));
    }

    /*printf("A(%d,%d) = %lg; y[%d] = %lg\n",j,j,Ajj,i,y[i]);*/

    for (ip = (L->colptr)[j] + 1; ip < (L->colptr)[j+1]; ip++) {
      i = (L->rowind)[ip];
      Aij = (L->taucs_values)[ip];

      /*x[i] -= y[j]*Aij;*/
      x[i] = taucs_sub( x[i], taucs_mul( y[j],Aij ));
    }
  }
  
  /* Solve D y = y */

  for (j=0; j<n; j++) {

    /* we put diagonal elements first */
    ip = (L->colptr)[j];
    i = (L->rowind)[ip];
    assert (i==j);
    Ajj = (L->taucs_values)[ip];
    
    /* y[j] = y[j] / Ajj; */
    y[j] = taucs_div( y[j] , Ajj );
  }

  /* Solve L^T x = y */

  /* taucs_printf("taucs_ccs_solve_ldlt: solving L^T x = y\n");*/

  for (i=n-1; i>=0; i--) {

    for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++) {
      j = (L->rowind)[jp];
      /* Aij = (L->taucs_values)[jp]; */
      Aij = taucs_conj( (L->taucs_values)[jp] );
      /*y[i] -= x[j]*Aij;*/
      y[i] = taucs_sub( y[i] , taucs_mul( x[j],Aij ));
    }

#if 0
    jp = (L->colptr)[i];
    j = (L->rowind)[jp];
    /*Aii = 1.0;*/
    Aii = taucs_one;

    /* x[i] = y[i] / Aii;*/
    x[i] = taucs_div( y[i] , Aii );
#else
    x[i] = y[i];
#endif

    if (taucs_isnan(x[i]) || taucs_isinf(x[i]))
	      taucs_printf("symccs_solve_ldlt: inf/nan in row %d (LT)\n",i);

    /*printf("A(%d,%d) = %lg; x[%d] = %lg\n",i,i,Aii,i,x[i]); */
  }

  taucs_free(y);

  return 0;
}

#endif /*#ifndef TAUCS_CORE_GENERAL*/

#ifdef TAUCS_CORE_GENERAL
int 
taucs_ccs_solve_schur(taucs_ccs_matrix* L,
		      taucs_ccs_matrix* schur_comp,
		      int    (*schur_precond_fn)(void*,void* x,void* b),
		      void*  schur_precond_args,
		      int    maxits,
		      double convratio,
		      void* x, void* b)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (L->flags & TAUCS_DOUBLE)
    return taucs_dccs_solve_schur(L,
				  schur_comp,
				  schur_precond_fn,
				  schur_precond_args,
				  maxits,
				  convratio,
				  (taucs_double*)x,(taucs_double*)b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (L->flags & TAUCS_SINGLE)
    return taucs_sccs_solve_schur(L,
				  schur_comp,
				  schur_precond_fn,
				  schur_precond_args,
				  maxits,convratio,
				  (taucs_single*)x,(taucs_single*)b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_solve_schur(L,
				  schur_comp,
				  schur_precond_fn,
				  schur_precond_args,
				  maxits,convratio,
				  (taucs_dcomplex*)x,
				  (taucs_dcomplex*)b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_solve_schur(L,
				  schur_comp,
				  schur_precond_fn,
				  schur_precond_args,
				  maxits,convratio,
				  (taucs_scomplex*)x,(taucs_scomplex*)b);
#endif
  
  assert(0);
  return -1;
}

int
taucs_ccs_solve_llt(void* vL, void* x, void* b)
{
  taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (L->flags & TAUCS_DOUBLE)
    return taucs_dccs_solve_llt(L,(taucs_double*) x, (taucs_double*) b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (L->flags & TAUCS_SINGLE)
    return taucs_sccs_solve_llt(L,(taucs_single*) x, (taucs_single*) b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_solve_llt(L,(taucs_dcomplex*) x, (taucs_dcomplex*) b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_solve_llt(L,(taucs_scomplex*) x, (taucs_scomplex*) b);
#endif
  
  assert(0);
  return -1;
}

int
taucs_ccs_solve_ldlt(void* vL, void* x, void* b)
{
  taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (L->flags & TAUCS_DOUBLE)
    return taucs_dccs_solve_ldlt(L,(taucs_double*) x, (taucs_double*) b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (L->flags & TAUCS_SINGLE)
    return taucs_sccs_solve_ldlt(L,(taucs_single*) x, (taucs_single*) b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_solve_ldlt(L,(taucs_dcomplex*) x, (taucs_dcomplex*) b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (L->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_solve_ldlt(L,(taucs_scomplex*) x, (taucs_scomplex*) b);
#endif

  /*omer*/
  assert(0);
  return -1;

}
#endif /*TAUCS_CORE_GENERAL*/


/*********************************************************/
/*                                                       */
/*********************************************************/
