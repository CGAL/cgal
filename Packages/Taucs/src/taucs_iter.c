/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "taucs.h"

#ifdef TAUCS_CORE_DOUBLE

/*********************************************************/
/* utilities                                             */
/*********************************************************/
/*extern int _isnan(double);*/

static double dotprod(int n, double* v, double* u)
{
  double x;
  int i;

  for (i=0, x=0.0; i<n; i++) x += v[i]*u[i];

  return x;
}

static double twonorm(int n, double* v)
{
  /*
  double norm;
  int i;

  for (i=0, norm=0.0; i<n; i++) norm += v[i]*v[i];

  norm = sqrt(norm);
  return norm;
  */

  double ssq, scale, absvi;/*norm omer*/
  int i;

  if (n==1) return fabs(v[0]);

  scale = 0.0;
  ssq   = 1.0;

  for (i=0; i<n; i++) {
    if ( v[i] != 0 ) {
      absvi = fabs(v[i]);
      if (scale < absvi) {
	ssq   = 1.0 + ssq * (scale/absvi)*(scale/absvi);
	scale = absvi;
      } else
	ssq   = ssq + (absvi/scale)*(absvi/scale);
    }
  }
  return scale * sqrt( ssq );
}

/*********************************************************/
/* conjugate gradients                                   */
/*********************************************************/

int 
taucs_conjugate_gradients(taucs_ccs_matrix* A,
			  int               (*precond_fn)(void*,void* x,void* b),
			  void*             precond_args,
			  void*             vX,
			  void*             vB,
			  int               itermax,
			  double            convergetol
			  )
{
  double* X = (double*) vX;
  double* B = (double*) vB;
  double *P, *R, *Q, *Z ;
  double Alpha, Beta, Rho, Init_norm, ratio, Res_norm, Rtmp ;
  double Rho0 = 0.0; /* warning */
  /*double t1, t2,  cpus[9] ; omer*/
  /*
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0} ;
  */
  /*  double Tiny = 0.0;*/
  double Tiny = 0.1e-28;
  int    Iter;
  /*int    stats[6] ; omer*/
  int    i,n;

#define RESVEC_NO
#ifdef RESVEC
  FILE* f;
  double* resvec = (double*) taucs_malloc((itermax+2) * sizeof(double));
  assert(resvec);
  for (i=0; i<=itermax; i++) {
    /*double inf = 1.0/0.0; omer*/
    double nan = taucs_get_nan()/*inf - inf; omer*/
    assert(taucs_isnan(nan));
    resvec[i] = nan;
  }
#endif

  n = A->n;
 
  P = (double*) taucs_malloc(n * sizeof(double));
  R = (double*) taucs_malloc(n * sizeof(double));
  Q = (double*) taucs_malloc(n * sizeof(double));
  Z = (double*) taucs_malloc(n * sizeof(double));

#define TAUCS_REMOVE_CONST_NO
#ifdef TAUCS_REMOVE_CONST
    {
      double s;
      for (i=0, s=0.0; i<n; i++) s += B[i];
      for (i=0, s=0.0; i<n; i++) B[i] -= s;
    }
#endif

  /*
  for (i=0; i<n; i++) X[i] = 0;
  for (i=0; i<n; i++) R[i] = B[i];
  */

  taucs_ccs_times_vec(A,X,R);
  for (i=0; i<n; i++) R[i] = B[i] - R[i];

  Res_norm = Init_norm = twonorm(n,R);
  printf("two norm of initial residual %.2e\n",Init_norm);
  if ( Init_norm == 0.0 ) Init_norm = 1.0;
  ratio = 1.0;
 
  Iter = 0;
 
#ifdef RESVEC
  resvec[Iter] = Res_norm;
#endif

  while ( ratio > convergetol && Iter <= itermax ) {
    Iter++;
    
    if (precond_fn)
      (*precond_fn)(precond_args,Z,R);
    else
      for (i=0; i<n; i++) Z[i] = R[i];

    for (i=0,Rho=0.0; i<n; i++) Rho += R[i] * Z[i];

    if ( Iter == 1 ) {
      for (i=0; i<n; i++) P[i] = Z[i];
    } else {
      Beta = Rho /(Rho0 + Tiny);
      for (i=0; i<n; i++) P[i] = Z[i] + Beta * P[i];
    };
 
    taucs_ccs_times_vec(A,P,Q); /* Q = A*P */

    for (i=0,Rtmp=0.0; i<n; i++) Rtmp += P[i] * Q[i];

    Alpha = Rho/(Rtmp+Tiny);

    for (i=0; i<n; i++) X[i] = X[i] + Alpha * P[i];

    for (i=0; i<n; i++) R[i] = R[i] - Alpha * Q[i];

#ifdef TAUCS_REMOVE_CONST
    {
      double s;
      for (i=0, s=0.0; i<n; i++) s += R[i];
      for (i=0, s=0.0; i<n; i++) R[i] -= s;
    }
#endif


    Rho0  = Rho;

    Res_norm = twonorm(n,R);

#if 0
    taucs_ccs_times_vec(A,X,R);
    for (i=0; i<n; i++) R[i] -= B[i];
    Res_norm = twonorm(n,R);
#endif

#ifdef RESVEC
  resvec[Iter] = Res_norm;
#endif

    ratio = Res_norm/Init_norm;
    if (Iter % 25 == 0) 
      taucs_printf("cg: n=%d at iteration %d the convergence ratio is %.2e, Rnorm %.2e\n", 
		   A->n,Iter, ratio,Res_norm) ;
  }
  if (Iter > 0) {
    taucs_printf("cg: n=%d iterations = %d Reduction in residual norm %.2e, Rnorm %.2e\n", 
		 A->n,Iter,ratio,Res_norm) ;
    taucs_ccs_times_vec(A,X,R);
    for (i=0; i<n; i++) R[i] = B[i] - R[i];
    taucs_printf("cg: true residual norm %.2e\n",twonorm(n,R));
  }

  taucs_free(P) ;
  taucs_free(R) ;
  taucs_free(Q) ;
  taucs_free(Z) ;
 
#ifdef RESVEC
  f=fopen("resvec","a");
  assert(f);
  for (i=0; i<=itermax && !taucs_isnan(resvec[i]); i++) {
    fprintf(f,"%.3e\n",resvec[i]);
  }
  fclose(f);
  taucs_free(resvec);
#endif

  return 0; 
}                                                                             

/*********************************************************/
/* minres                                                */
/*********************************************************/

int 
taucs_minres(taucs_ccs_matrix*  A,
	     int                (*precond_fn)(void*,void* x,void* b),
	     void*              precond_args,
	     void*              vX,
	     void*              vB,
	     int                itermax,
	     double             convergetol)
{
  double* X = (double*) vX;
  double* B = (double*) vB;

  double *Xcg, *R, *V, *VV, *Vold, *Volder, *M, *Mold, *Molder;
  double tolb, normr, alpha, beta, beta1, betaold;
  double gamma, gammabar, delta, deltabar, epsilon;
  double cs,sn,snprod, numer, denom;
  int    Iter;
  int    i,n;

  n = A->n;
 
  R      = (double*) taucs_malloc(n * sizeof(double));
  Xcg    = (double*) taucs_malloc(n * sizeof(double));
  VV     = (double*) taucs_malloc(n * sizeof(double));
  V      = (double*) taucs_malloc(n * sizeof(double));
  Vold   = (double*) taucs_malloc(n * sizeof(double));
  Volder = (double*) taucs_malloc(n * sizeof(double));
  M      = (double*) taucs_malloc(n * sizeof(double));
  Mold   = (double*) taucs_malloc(n * sizeof(double));
  Molder = (double*) taucs_malloc(n * sizeof(double));

  tolb = convergetol * twonorm(n,B);
  taucs_printf("minres: residual convergence tolerance %.1e\n",tolb);
 
  for (i=0; i<n; i++) X[i] = 0;    /* x = 0 */
  for (i=0; i<n; i++) R[i] = B[i]; /* r = b-A*x */

  normr = twonorm(n,R);
  if ( normr == 0.0 ) {
    taucs_printf("minres: initial residual == 0\n");
    return -1;
  }

  for (i=0; i<n; i++) V[i]    = R[i];    /* v = r */
  for (i=0; i<n; i++) Vold[i] = R[i];    /* vold = r */
  
  if (precond_fn)
    (*precond_fn)(precond_args,V,Vold);
  else
    for (i=0; i<n; i++) V[i] = Vold[i];
  
  beta1 = dotprod(n,Vold,V);
  if (beta1 < 0.0) {
    taucs_printf("minres: error (1)\n");
    return -1;
  }
  beta1 = sqrt(beta1);

  { int flag = 0;
    for (i=0; i<n; i++) {
      if (taucs_isnan(V[i]) && flag < 10) 
	taucs_printf("minres: V has nan's in position %d\n",i);
      flag++;
    }
  }


  snprod = beta1;
  taucs_printf(">>> %e %e %e\n",beta1,snprod,normr);


  for (i=0; i<n; i++) VV[i] = V[i] / beta1;
  
  taucs_ccs_times_vec(A,VV,V); /* V = A*VV */
  
  alpha = dotprod(n,VV,V);
  
  for (i=0; i<n; i++) V[i] -= (alpha/beta1) * Vold[i];
  
  /* local reorthogonalization */

  numer = dotprod(n,VV,V);
  denom = dotprod(n,VV,VV);

  for (i=0; i<n; i++) V[i] -= (numer/denom) * VV[i];

  for (i=0; i<n; i++) Volder[i] = Vold[i];
  for (i=0; i<n; i++) Vold[i]   = V[i];
  
  if (precond_fn)
    (*precond_fn)(precond_args,V,Vold);
  else
    for (i=0; i<n; i++) V[i] = Vold[i];
  
  betaold = beta1;
  beta = dotprod(n,Vold,V);
  if (beta < 0.0) {
    taucs_printf("minres: error (2)\n");
    return -1;
  }
  beta = sqrt(beta);
  
  gammabar = alpha;
  epsilon = 0.0;
  deltabar = beta;
  gamma = sqrt(gammabar*gammabar + beta*beta);


  for (i=0; i<n; i++) Mold[i] = 0.0;
  for (i=0; i<n; i++) M[i]    = VV[i] / gamma;

  cs = gammabar / gamma;
  sn = beta / gamma;


  for (i=0; i<n; i++) X[i] += snprod*cs*M[i];
  snprod = snprod * sn;

  /* generate CG iterates */
  for (i=0; i<n; i++) Xcg[i] = X[i] + snprod*(sn/cs)*M[i];

  /* compute residual again */
  
  taucs_ccs_times_vec(A,X,R); 
  for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
  normr = twonorm(n,R);

  taucs_printf("minres: starting iterations, residual norm is %.1e\n",normr);
  
  for ( Iter=1; Iter <= itermax; Iter++ ) {

    for (i=0; i<n; i++) VV[i] = V[i] / beta;
    taucs_ccs_times_vec(A,VV,V); 
    for (i=0; i<n; i++) V[i] -= (beta/betaold) * Volder[i];
    alpha = dotprod(n,VV,V);
    for (i=0; i<n; i++) V[i] -= (alpha/beta) * Vold[i];

    for (i=0; i<n; i++) Volder[i] = Vold[i];
    for (i=0; i<n; i++) Vold  [i] = V   [i];
    
    if (precond_fn)
      (*precond_fn)(precond_args,V,Vold);
    else
      for (i=0; i<n; i++) V[i] = Vold[i];

    betaold = beta;
    beta = dotprod(n,Vold,V);
    if (beta < 0.0) {
      taucs_printf("minres: error (3)\n");
      return -1;
    }
    beta = sqrt(beta);

    delta = cs*deltabar + sn*alpha;
    for (i=0; i<n; i++) Molder[i] = Mold[i];
    for (i=0; i<n; i++) Mold  [i] = M   [i];
    for (i=0; i<n; i++) M[i] = VV[i] - delta*Mold[i] - epsilon*Molder[i];
    gammabar = sn*deltabar - cs*alpha;
    epsilon = sn*beta;
    deltabar = -cs*beta;
    gamma = sqrt(gammabar*gammabar + beta*beta);
    for (i=0; i<n; i++) M[i] = M[i]/ gamma;
    cs = gammabar / gamma;
    sn = beta / gamma;

    /* stagnation test; skipped */
    
    for (i=0; i<n; i++) X[i] += snprod*cs*M[i];
    snprod = snprod*sn;
    for (i=0; i<n; i++) Xcg[i] = X[i] + snprod*(sn/cs)*M[i];
    
    if (precond_fn) {
      taucs_ccs_times_vec(A,X,R); 
      for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
      normr = twonorm(n,R);
    } else {
      normr = fabs(snprod); 
      if (normr <= tolb) {
	/* double check */
	taucs_ccs_times_vec(A,X,R); 
	for (i=0; i<n; i++) R[i] = B[i] - R[i];  /* r = b - A*x */
	normr = twonorm(n,R);
      }
    }

    if (Iter > -1)
      taucs_printf("minres: n=%d iterations = %d residual norm %12.4e\n", A->n,Iter,normr);

    if (normr <= tolb) break;
  }

  taucs_printf("minres: done. n=%d iterations = %d residual norm %12.4e\n", A->n,Iter,normr);
 
  taucs_free(Molder) ;
  taucs_free(Mold) ;
  taucs_free(M) ;
  taucs_free(Volder) ;
  taucs_free(Vold) ;
  taucs_free(V) ;
  taucs_free(VV) ;
  taucs_free(Xcg) ;
  taucs_free(R) ;
 
  return 0; 
}                                                                             

#endif /* TAUCS_CORE_DOUBLE */

