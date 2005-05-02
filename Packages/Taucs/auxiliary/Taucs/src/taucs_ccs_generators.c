/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "taucs.h"

#define RNDM ((double)rand()/(double)RAND_MAX);

/*ifndef added omer*/
#ifndef max
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif
#ifndef mod
#define mod(x,n) ((x) % (n))
#endif
/*********************************************************/
/*                                                       */
/*********************************************************/

#ifdef TAUCS_CORE_DOUBLE

taucs_ccs_matrix* taucs_ccs_generate_mesh2d_negative(int n)
{
  taucs_ccs_matrix* m;
  int         N;
  int         nnz;
  int         x,y,i,j,ip;

  taucs_printf("generate_mesh2d_negative: starting\n");

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("generate_mesh2d_negative: out of memory (1)\n");
    return NULL; 
  }

  N   = n*n;
  nnz = 4*N;

  m->n      = N;
  m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("generate_mesh2d_negative: out of memory (4): ncols=%d nnz=%d\n",N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/);
    return NULL; 
  }

  ip = 0;
  for (y=0; y<n; y++) {
    for (x=0; x<n; x++) {
      j = x + y*n; /* convert mesh (x,y) location to index in vector */
      /*printf("column %d xy %d,%d starts at %d\n",j,x,y,ip);*/
      (m->colptr)[j] = ip;

      i=mod(x+1,n) + (y  )*n      ; if (i>j) { (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }
      i=(x  )      + mod(y+1,n)*n ; if (i>j) { (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=+100.0; ip++; }
      i=mod(x+n-1,n) + (y  )*n      ; if (i>j) { (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }
      i=(x  )      + mod(y+n-1,n)*n ; if (i>j) { (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=+100.0; ip++; }

      /* i=mod(x+1,n) + mod(y+1,n)*n ; (m->rowind)[ip]=i; (m->taucs_values)[ip]=+1.0; ip++;  */
      /* i=mod(x+2,n) + (y  )*n      ; (m->rowind)[ip]=i; (m->taucs_values)[ip]=.0625; ip++;  */
      /* i=(x  )      + mod(y+2,n)*n ; (m->rowind)[ip]=i; (m->taucs_values)[ip]=.0625; ip++;  */

      i=(x  )+(y  )*n; (m->rowind)[ip]=i;
      /* (m->taucs_values)[ip]= 4.25; if (x==0 && y==0) (m->taucs_values)[ip] += 1; to make it nonsingular  */
      (m->values.d/*taucs_values*/)[ip]= 202.0; if (x==0 && y==0) (m->values.d/*taucs_values*/)[ip] += 1; /* to make it nonsingular */ 
      ip++; 
    }
  }

  (m->colptr)[N] = ip;

  taucs_printf("generate_mesh2d_negative: done: ncols=%d nnz=%d\n",N,ip);

  return m;
}


taucs_ccs_matrix* 
taucs_ccs_generate_mesh2d(int n,char *which)
{
  taucs_ccs_matrix* m;
  int         N;
  int         nnz;
  int         x,y,i,j,ip;
  double jump = 100;

  taucs_printf("taucs_ccs_generate_mesh2d: starting\n");

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("generate_mesh2d: out of memory (1)\n");
    return NULL; 
  }

  N   = n*n;
  nnz = 3*N;

  m->n      = N;
  m->m      = N;
  m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("taucs_ccs_generate_mesh2d: out of memory: ncols=%d nnz=%d\n",N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/);
    return NULL; 
  }

  ip = 0;
  for (y=0; y<n; y++) {
    for (x=0; x<n; x++) {
      j = x + y*n; /* convert mesh (x,y) location to index in vector */
      /*printf("column %d xy %d,%d starts at %d\n",j,x,y,ip);*/
      (m->colptr)[j] = ip;
      /* if (x < n-1) { i=(x+1)+(y  )*n; (m->rowind)[ip]=i; (m->taucs_values)[ip]=-1.0; ip++; } */
      
      if (!strcmp(which,"anisotropic_y")) {
	if (y < n-1) { i=(x  )+(y+1)*n; (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-jump; ip++; }
      } else 
	if (y < n-1) { i=(x  )+(y+1)*n; (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }

      if (!strcmp(which,"anisotropic_x")) {
	if (x < n-1) { i=(x+1)+(y  )*n; (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-jump; ip++; }
      } else 
	if (x < n-1) { i=(x+1)+(y  )*n; (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }

      if (!strcmp(which,"anisotropic_y")) 
	{ 
	  i=(x  )+(y  )*n; (m->rowind)[ip]=i; 
	  (m->values.d/*taucs_values*/)[ip]= 0.0; 
	  if (x > 0)   (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (y > 0)   (m->values.d/*taucs_values*/)[ip] += jump; 
	  if (x < n-1) (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (y < n-1) (m->values.d/*taucs_values*/)[ip] += jump; 
	  if (x==0 && y==0) (m->values.d/*taucs_values*/)[ip] += 1.0; /* to make it nonsingular */
	  ip++; 
	}
      else if (!strcmp(which,"anisotropic_x")) 
	{ 
	  i=(x  )+(y  )*n; (m->rowind)[ip]=i; 
	  (m->values.d/*taucs_values*/)[ip]= 0.0; 
	  if (x > 0)   (m->values.d/*taucs_values*/)[ip] += jump; 
	  if (y > 0)   (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (x < n-1) (m->values.d/*taucs_values*/)[ip] += jump; 
	  if (y < n-1) (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (x==0 && y==0) (m->values.d/*taucs_values*/)[ip] += 1.0; /* to make it nonsingular */
	  ip++; 
	}
      else if (!strcmp(which,"dirichlet"))
	{
	  i=(x  )+(y  )*n; (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]= 4.0; 
	  ip++; 
	}
      else /* neumann */
	{ 
	  i=(x  )+(y  )*n; (m->rowind)[ip]=i; 
	  (m->values.d/*taucs_values*/)[ip]= 0.0; 
	  if (x > 0)   (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (y > 0)   (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (x < n-1) (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (y < n-1) (m->values.d/*taucs_values*/)[ip] += 1.0; 
	  if (x==0 && y==0) (m->values.d/*taucs_values*/)[ip] += 1.0; /* to make it nonsingular */
	  ip++; 
	}
	
    }
  }
  (m->colptr)[N] = ip;

  taucs_printf("taucs_ccs_generate_mesh2d: done, ncols=%d nnz=%d\n",N,ip);

  /*
  for (j=0; j<N; j++) {
    for (ip=(m->colptr)[j]; ip < (m->colptr)[j+1]; ip++) {
      i = (m->rowind)[ip];
      taucs_printf("<%d %d %lg>\n",i,j,m->taucs_values[ip]);
    }
  }
  */

  return m;
}

taucs_ccs_matrix* 
taucs_ccs_generate_mesh3d(int X, int Y, int Z)
{
  taucs_ccs_matrix* m;
  int         N;
  int         nnz;
  int         x,y,z,i,j,ip;

  taucs_printf("taucs_ccs_generate_mesh3d: starting\n");

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_generate_mesh3d: out of memory\n");
    return NULL; 
  }

  N   = X*Y*Z;
  nnz = 4*N;

  m->n      = N;
  m->m      = N;
  m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  /*m->indshift = 0;*/
  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("taucs_ccs_generate_mesh3d: out of memory: ncols=%d nnz=%d\n",N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/);
    return NULL; 
  }

  ip = 0;
  for (z=0; z<Z; z++) {
    for (y=0; y<Y; y++) {
      for (x=0; x<X; x++) {
	j = z*X*Y + y*X + x; 
	/*printf("column %d xy %d,%d starts at %d\n",j,x,y,ip);*/
	(m->colptr)[j] = ip;
	if (x < X-1) { i=(z  )*X*Y+(y  )*X+(x+1); (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }
	if (y < Y-1) { i=(z  )*X*Y+(y+1)*X+(x  ); (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }
	if (z < Z-1) { i=(z+1)*X*Y+(y  )*X+(x  ); (m->rowind)[ip]=i; (m->values.d/*taucs_values*/)[ip]=-1.0; ip++; }
	             { 
		       i=(z  )*X*Y+(y  )*X+(x  ); (m->rowind)[ip]=i; 
		       (m->values.d/*taucs_values*/)[ip]= 0.0; 
		       if (x < X-1) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (y < Y-1) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (z < Z-1) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (x > 0  ) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (y > 0  ) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (z > 0  ) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       if (x==0 && y==0 && z==0) (m->values.d/*taucs_values*/)[ip] += 1.0;
		       ip++; 
		     }
	/* { i=(z  )*X*Y+(y  )*X+(x  ); (m->rowind)[ip]=i; (m->taucs_values)[ip]= 6.0; ip++; } */
      }
    }
  }
  (m->colptr)[N] = ip;

  taucs_printf("taucs_ccs_generate_mesh3d: done, ncols=%d nnz=%d\n",N,ip);

  return m;
}

taucs_ccs_matrix* 
taucs_ccs_generate_dense(int M, int N, int flags)
{
  taucs_ccs_matrix* m;
  int         nnz;
  int         i,j,ip;/* x,y omer*/

  taucs_printf("taucs_ccs_generate_dense: starting\n");

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_generate_dense: out of memory\n");
    return NULL; 
  }

  m->m      = N;
  m->n      = N;
  if (flags & TAUCS_SYMMETRIC) {
    nnz = N*(N+1)/2;
    m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  } else {
    nnz = N*N;
    m->flags  =  TAUCS_DOUBLE;
  }

  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("taucs_ccs_generate_dense: out of memory: nrows=%d ncols=%d nnz=%d\n",M,N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/);
    return NULL; 
  }

  ip = 0;
  for (j=0; j<N; j++) {
    (m->colptr)[j] = ip;
    if (flags & TAUCS_SYMMETRIC) {
      for (i=j; i<N; i++) {
	(m->rowind)[ip]=i; 
	(m->values.d/*taucs_values*/)[ip]=RNDM; 
	ip++;
      }
    } else {
      for (i=0; i<M; i++) {
	(m->rowind)[ip]=i; 
	(m->values.d/*taucs_values*/)[ip]=RNDM; 
	ip++;
      }
    }
  }
  (m->colptr)[N] = ip;

  taucs_printf("taucs_ccs_generate_dense: done, nrows=%d ncols=%d nnz=%d\n",M,N,ip);

  return m;
}

/* random resistor networks */

int recursive_visit(int i, 
		    int* neighbors[], 
		    int degree[], 
		    int visited[])
{
  int j,jp,count;
  visited[i] = 1;
  count = 1;
  for (jp=0; jp<degree[i]; jp++) {
    j = neighbors[i][jp];
    if (! visited[j] ) count += recursive_visit(j,neighbors,degree,visited);
  }
  return count;
}

taucs_ccs_matrix* 
taucs_ccs_generate_rrn(int X, int Y, int Z, double drop_probability, double rmin)
{
  taucs_ccs_matrix* m;
  taucs_ccs_matrix* l;
  int         N;
  int         nnz;
  int         x,y,z,i,j,k,ip,jp;
  double*     D; /* contributions to future diagonal elements */

  int**       neighbors;
  int*        degree;
  int*        visited;
  int*        reps;
  int         ncomponents;

  int         largest;
  int         largest_rep = 0; /* warning */

  taucs_printf("taucs_ccs_generate_rrn: starting (%d %d %d %.4e %.4e)\n",
	       X,Y,Z,drop_probability,rmin);

  if (drop_probability > 1.0 || drop_probability < 0.0) {
    taucs_printf("taucs_ccs_generate_rrn: drop probability (%lg) must be in [0,1], setting to 0\n",
		 drop_probability);
    drop_probability = 0.0;
  }

  if (rmin > 1.0 || rmin <= 0.0) {
    taucs_printf("taucs_ccs_generate_rrn: rmin (%lg) must be in (0,1], setting to 1\n",
		 rmin);
    rmin = 1.0;
  }

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_generate_rrn: out of memory\n");
    return NULL; 
  }

  N   = X*Y*Z;
  nnz = 4*N;   /* this is an upper bound */

  m->n      = N;
  m->m      = N;
  m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  /*m->indshift = 0;*/
  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  D         = (double*) taucs_malloc(N         * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind) || !D) {
    taucs_printf("taucs_ccs_generate_rrn: out of memory: ncols=%d nnz=%d\n",N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/); taucs_free(D);
    return NULL; 
  }

  for (i=0; i<N; i++) D[i] = 0.0;

  ip = 0;
  for (z=0; z<Z; z++) {
    for (y=0; y<Y; y++) {
      for (x=0; x<X; x++) {
	int j, je, jw, js, jn, ju, jd; /* indices for up, down, east, west, south, north */
	int jp; /* pointer to the diagonal value */
	double v;

	j  = z*X*Y + y*X + x;
	jw = (x > 0  ) ? (z  )*X*Y + (y  )*X + (x-1) : (z  )*X*Y + (y  )*X + (X-1) ;
	je = (x < X-1) ? (z  )*X*Y + (y  )*X + (x+1) : (z  )*X*Y + (y  )*X + (0  ) ;
	js = (y > 0  ) ? (z  )*X*Y + (y-1)*X + (x  ) : (z  )*X*Y + (Y-1)*X + (x  ) ;
	jn = (y < Y-1) ? (z  )*X*Y + (y+1)*X + (x  ) : (z  )*X*Y + (0  )*X + (x  ) ;
	jd = (z > 0  ) ? (z-1)*X*Y + (y  )*X + (x  ) : (Z-1)*X*Y + (y  )*X + (x  ) ;
	ju = (z < Z-1) ? (z+1)*X*Y + (y  )*X + (x  ) : (0  )*X*Y + (y  )*X + (x  ) ;

	jw = (x > 0  ) ? (z  )*X*Y + (y  )*X + (x-1) : -1;
	je = (x < X-1) ? (z  )*X*Y + (y  )*X + (x+1) : -1;
	js = (y > 0  ) ? (z  )*X*Y + (y-1)*X + (x  ) : -1;
	jn = (y < Y-1) ? (z  )*X*Y + (y+1)*X + (x  ) : -1;
	jd = (z > 0  ) ? (z-1)*X*Y + (y  )*X + (x  ) : -1;
	ju = (z < Z-1) ? (z+1)*X*Y + (y  )*X + (x  ) : -1;

	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) jw = -1;
	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) je = -1;
	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) js = -1;
	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) jn = -1;
	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) ju = -1;
	if ( ((double)rand() / (double)RAND_MAX) < drop_probability) jd = -1;

	/*printf("xyz=%d %d %d    j's=%d %d %d %d %d %d %d\n",x,y,z,j,jw,je,js,jn,jd,ju);*/
	/*printf("column %d xy %d,%d starts at %d\n",j,x,y,ip);*/
	(m->colptr)[j] = ip;
	jp = ip;
	
	/*printf("j=%d D[]=%lf\n",j,D[j]);*/

	(m->rowind)[ip]= j;
	/*
	if (x==0 && y==0 && z==0)
	  (m->taucs_values)[ip]= 1.0;
	else 
	*/
	(m->values.d/*taucs_values*/)[ip]= D[j]; 
	ip++;

	if (jw != j != -1) {
	  if (jw > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jw;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jw] -= v;
	  }
	}

	if (je != j && je != jw && je != -1) {
	  if (je > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = je;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[je] -= v;
	  }
	}

	if (js != j && js != -1) {
	  if (js > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = js;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[js] -= v;
	  }
	}

	if (jn != j && jn != js && jn != -1) {
	  if (jn > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jn;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jn] -= v;
	  }
	}

	if (ju != j && ju != -1) {
	  if (ju > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = ju;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[ju] -= v;
	  }
	}

	if (jd != j && jd != ju && jd != -1) {
	  if (jd > j) {
	    v = -1.0;
	    v = ((double)rand()/(double)RAND_MAX) > 0.99 ? -1.0 : -rmin;
	    v = -( rmin + (((double)rand()/(double)RAND_MAX) * (1.0-rmin)) );
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jd;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jd] -= v;
	  }
	}

      }
    }
  }
  taucs_free(D);
  (m->colptr)[N] = ip;

  taucs_printf("taucs_ccs_generate_rrn: done, ncols=%d allocated nnz=%d real nnz=%d\n",
	       N,nnz,ip);


  neighbors = (int**) taucs_malloc(N * sizeof(int*));
  degree  = (int*) taucs_malloc(N * sizeof(int));
  visited = (int*) taucs_malloc(N * sizeof(int));
  reps    = (int*) taucs_malloc(N * sizeof(int));

  for (i=0; i<N; i++) degree[i] = 0;

  for (j=0; j<N; j++) {
    for (ip=(m->colptr)[j]; ip<(m->colptr)[j+1]; ip++) {
      i = (m->rowind)[ ip ];
      if (i != j) {
	degree[i]++;
	degree[j]++;
      }
    }
  }


  for (i=0; i<N; i++) {
    neighbors[i] = (int*) taucs_malloc(degree[i] * sizeof(int));
    visited[i] = 0;
  }

  for (j=0; j<N; j++) {
    for (ip=(m->colptr)[j]; ip<(m->colptr)[j+1]; ip++) {
      i = (m->rowind)[ ip ];
      if (i != j) {
	neighbors[i][visited[i]] = j;
	neighbors[j][visited[j]] = i;
	assert(visited[i] < degree[i]);
	assert(visited[j] < degree[j]);
	visited[i]++;
	visited[j]++;
      }
    }
  }

  for (i=0; i<N; i++) visited[i] = 0;
  ncomponents = 0;
  largest = -1;
  for (i=0; i<N; i++) {
    if (visited[i] == 0) {
      int count;
      reps[ncomponents] = i;
      ncomponents++;
      count = recursive_visit(i,neighbors,degree,visited);
      if (count > largest) {
	largest = count;
	largest_rep = i;
      }
      /*printf("new connected component vertex %d, size=%d\n",i,count);*/
    }
  }
  for (i=0; i<ncomponents; i++) {
    j = reps[i];
    /*printf("rep[%d] = %d\n",i,j);*/
    (m->values.d/*taucs_values*/)[ (m->colptr)[j] ] += 1.0;
  }
  printf("found %d components, largest is %d, rep is %d\n",ncomponents,largest,largest_rep);
  printf("found %d components\n",ncomponents);

  for (i=0; i<N; i++) visited[i] = 0;
  (void) recursive_visit(largest_rep,neighbors,degree,visited);
  
  /* we now reuse the degree and reps vectors */

  for (i=0; i<N; i++) degree[i] = reps[i] = -1;
  j = 0;
  for (i=0; i<N; i++) {
    if (visited[i]) {
      degree[i] = j;
      reps[j] = i;
      j++;
    }
  }

  l = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!l) { 
    taucs_printf("taucs_ccs_generate_rrn: out of memory\n");
    return NULL; 
  }

  nnz = (m->colptr)[N];   /* this is an upper bound */

  l->n      = largest;
  l->m      = largest;
  l->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  l->colptr = (int*)    taucs_malloc((largest+1) * sizeof(int));
  l->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  l->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  k = 0;
  for (jp=0; jp<N; jp++) {
    int iip;
    j = degree[jp];
    if (j == -1) continue;
    assert(j < largest);
    (l->colptr)[j] = k;
    for (iip=(m->colptr)[jp]; iip<(m->colptr)[jp+1]; iip++) {
      double v;
      ip = (m->rowind)[iip];
      v  = (m->values.d/*taucs_values*/)[iip];
      i = degree[ip];
      assert(i >= j);
      (l->rowind)[k] = i;
      (l->values.d/*taucs_values*/)[k] = v;
      k++;
    }
  }
  (l->colptr)[largest] = k;

  for (i=0; i<N; i++) taucs_free(neighbors[i]);
  taucs_free(visited);
  taucs_free(reps);
  taucs_free(degree);
  taucs_free(neighbors);

  taucs_ccs_free(m);

  return l;
}

double* taucs_vec_generate_continuous(int X, int Y, int Z, char* which)
{
  int x,y,z,j;/* i,k omer*/
  double* V;
  double dx,dy,dz;

  V = (double*) taucs_malloc( X*Y*Z * sizeof(double));
  if (!V) {
    taucs_printf("taucs_vec_generate_continuous: out of memory\n");
    return V;
  }

  for (z=0; z<Z; z++) {
    for (y=0; y<Y; y++) {
      for (x=0; x<X; x++) {
	double v;

	j  = z*X*Y + y*X + x;
	
	dx = (double) (x+1) / (double) X;
	dy = (double) (y+1) / (double) Y;
	dz = (double) (z+1) / (double) Z;

	v = (dx*dy*dz*(1.0-dx)*(1.0-dy)*(1.0-dz));
	v = v*v;
	v = v*exp(dx*dx*dy*dz);

	V[j] = v;
      }
    }
  }

  return V;
}

taucs_ccs_matrix* 
taucs_ccs_generate_discontinuous(int X, int Y, int Z, double jump)
{
  taucs_ccs_matrix* m;
  /*taucs_ccs_matrix* l; omer*/
  int         N;
  int         nnz;
  int         x,y,z,i,ip;/*j,k,jp omer*/
  double*     D; /* contributions to future diagonal elements */

  taucs_printf("taucs_ccs_generate_discontinuous: starting (%d %d %d %e)\n",
	       X,Y,Z,jump);


  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_generate_discontinuous: out of memory\n");
    return NULL; 
  }

  N   = X*Y*Z;
  nnz = 4*N;   /* this is an upper bound */

  m->n      = N;
  m->m      = N;
  m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  /*m->indshift = 0;*/
  m->colptr = (int*)    taucs_malloc((N+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->values.d/*taucs_values*/ = (double*) taucs_malloc(nnz       * sizeof(double));

  D         = (double*) taucs_malloc(N         * sizeof(double));

  if (!(m->colptr) || !(m->rowind) || !(m->rowind) || !D) {
    taucs_printf("taucs_ccs_generate_discontinuous: out of memory: ncols=%d nnz=%d\n",N,nnz);
    taucs_free(m->colptr); taucs_free(m->rowind); taucs_free(m->values.d/*taucs_values*/); taucs_free(D);
    return NULL; 
  }

  for (i=0; i<N; i++) D[i] = 0.0;

  ip = 0;
  for (z=0; z<Z; z++) {
    for (y=0; y<Y; y++) {
      for (x=0; x<X; x++) {
	int j, je, jw, js, jn, ju, jd; /* indices for up, down, east, west, south, north */
	int jp; /* pointer to the diagonal value */
	double v;
	int cj, cjw, cje, cjs, cjn, cjd, cju; /* which region? */

	j  = z*X*Y + y*X + x;
	jw = (x > 0  ) ? (z  )*X*Y + (y  )*X + (x-1) : (z  )*X*Y + (y  )*X + (X-1) ;
	je = (x < X-1) ? (z  )*X*Y + (y  )*X + (x+1) : (z  )*X*Y + (y  )*X + (0  ) ;
	js = (y > 0  ) ? (z  )*X*Y + (y-1)*X + (x  ) : (z  )*X*Y + (Y-1)*X + (x  ) ;
	jn = (y < Y-1) ? (z  )*X*Y + (y+1)*X + (x  ) : (z  )*X*Y + (0  )*X + (x  ) ;
	jd = (z > 0  ) ? (z-1)*X*Y + (y  )*X + (x  ) : (Z-1)*X*Y + (y  )*X + (x  ) ;
	ju = (z < Z-1) ? (z+1)*X*Y + (y  )*X + (x  ) : (0  )*X*Y + (y  )*X + (x  ) ;

	jw = (x > 0  ) ? (z  )*X*Y + (y  )*X + (x-1) : -1;
	je = (x < X-1) ? (z  )*X*Y + (y  )*X + (x+1) : -1;
	js = (y > 0  ) ? (z  )*X*Y + (y-1)*X + (x  ) : -1;
	jn = (y < Y-1) ? (z  )*X*Y + (y+1)*X + (x  ) : -1;
	jd = (z > 0  ) ? (z-1)*X*Y + (y  )*X + (x  ) : -1;
	ju = (z < Z-1) ? (z+1)*X*Y + (y  )*X + (x  ) : -1;

	/*printf("xyz=%d %d %d    j's=%d %d %d %d %d %d %d\n",x,y,z,j,jw,je,js,jn,jd,ju);*/
	/*printf("column %d xy %d,%d starts at %d\n",j,x,y,ip);*/
	(m->colptr)[j] = ip;
	jp = ip;
	
	/*printf("j=%d D[]=%lf\n",j,D[j]);*/

	(m->rowind)[ip]= j;
	/* Nonsingular Neumann */
	if (x==0 && y==0 && z==0)
	  (m->values.d/*taucs_values*/)[ip]= D[j] + 1.0;
	else 
	  (m->values.d/*taucs_values*/)[ip]= D[j]; 

	/* Singular Neumann */
	/*
	(m->taucs_values)[ip] = D[j];
	*/

	/* Dirichlet */
	/*
	(m->taucs_values)[ip] = D[j];
	if (x==0 || x==X-1) (m->taucs_values)[ip] += 1.0;
	if (y==0 || y==Y-1) (m->taucs_values)[ip] += 1.0;
	if (z==0 || z==Z-1) (m->taucs_values)[ip] += 1.0;
	*/

	ip++;

	cj  = ((x  ) >= X/8 && (x  ) < 7*X/8)
	   && ((y  ) >= Y/8 && (y  ) < 7*Y/8) 
	   && ((z  ) >= Z/8 && (z  ) < 7*Z/8);
	/*
	cj  = cj  && !(   ((x  ) >= 2*X/8 && (x  ) < 6*X/8)
		       && ((y  ) >= 2*Y/8 && (y  ) < 6*Y/8) 
		       && ((z  ) >= 2*Z/8 && (z  ) < 6*Z/8));
	*/
	cjw = ((x-1) >= X/8 && (x-1) < 7*X/8) 
	   && ((y  ) >= Y/8 && (y  ) < 7*Y/8) 
	   && ((z  ) >= Z/8 && (z  ) < 7*Z/8);
	/*
	cjw = cjw && !(   ((x-1) >= 2*X/8 && (x-1) < 6*X/8)
		       && ((y  ) >= 2*Y/8 && (y  ) < 6*Y/8) 
		       && ((z  ) >= 2*Z/8 && (z  ) < 6*Z/8));
	*/
	cje = ((x+1) >= X/8 && (x+1) < 7*X/8) 
	   && ((y  ) >= Y/8 && (y  ) < 7*Y/8) 
	   && ((z  ) >= Z/8 && (z  ) < 7*Z/8);
	/*
	cje = cje && !(   ((x+1) >= 2*X/8 && (x+1) < 6*X/8)
		       && ((y  ) >= 2*Y/8 && (y  ) < 6*Y/8) 
		       && ((z  ) >= 2*Z/8 && (z  ) < 6*Z/8));
	*/
	cjs = ((x  ) >= X/8 && (x  ) < 7*X/8) 
	   && ((y-1) >= Y/8 && (y-1) < 7*Y/8) 
	   && ((z  ) >= Z/8 && (z  ) < 7*Z/8);
	/*
	cjs = cjs && !(   ((x  ) >= 2*X/8 && (x  ) < 6*X/8)
		       && ((y-1) >= 2*Y/8 && (y-1) < 6*Y/8) 
		       && ((z  ) >= 2*Z/8 && (z  ) < 6*Z/8));
	*/
	cjn = ((x  ) >= X/8 && (x  ) < 7*X/8) 
	   && ((y+1) >= Y/8 && (y+1) < 7*Y/8) 
	   && ((z  ) >= Z/8 && (z  ) < 7*Z/8);
	/*
	cjn = cjn && !(   ((x  ) >= 2*X/8 && (x  ) < 6*X/8)
		       && ((y+1) >= 2*Y/8 && (y+1) < 6*Y/8) 
		       && ((z  ) >= 2*Z/8 && (z  ) < 6*Z/8));
	*/
	cjd = ((x  ) >= X/8 && (x  ) < 7*X/8) 
	   && ((y  ) >= Y/8 && (y  ) < 7*Y/8) 
	   && ((z-1) >= Z/8 && (z-1) < 7*Z/8);
	/*
	cjd = cjd && !(   ((x  ) >= 2*X/8 && (x  ) < 6*X/8)
		       && ((y  ) >= 2*Y/8 && (y  ) < 6*Y/8) 
		       && ((z-1) >= 2*Z/8 && (z-1) < 6*Z/8));
	*/
	cju = ((x  ) >= X/8 && (x  ) < 7*X/8) 
	   && ((y  ) >= Y/8 && (y  ) < 7*Y/8) 
	   && ((z+1) >= Z/8 && (z+1) < 7*Z/8);
	/*
	cju = cju && !(   ((x  ) >= 2*X/8 && (x  ) < 6*X/8)
		       && ((y  ) >= 2*Y/8 && (y  ) < 6*Y/8) 
		       && ((z+1) >= 2*Z/8 && (z+1) < 6*Z/8));
	*/

	if (jw != j && jw != -1) {
	  if (jw > j) {
	    v = -jump;
	    v = (x < X/8 || y < Y/8) ? -jump : -1.0;
	    v = (cj && cjw) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jw;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jw] -= v;
	  }
	}

	if (je != j && je != jw && je != -1) {
	  if (je > j) {
	    v = -jump;
	    v = ((x-1) < X/8 || y < Y/8) ? -jump : -1.0;
	    v = (cj && cje) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = je;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[je] -= v;
	  }
	}

	if (js != j && js != -1) {
	  if (js > j) {
	    v = -jump;
	    v = (y < Y/8 || x < X/8) ? -jump : -1.0;
	    v = (cj && cjs) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = js;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[js] -= v;
	  }
	}

	if (jn != j && jn != js && jn != -1) {
	  if (jn > j) {
	    v = -jump;
	    v = ((y-1) < Y/8 || x < X/8) ? -jump : -1.0;
	    v = (cj && cjn) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jn;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jn] -= v;
	  }
	}

	if (ju != j && ju != -1) {
	  if (ju > j) {
	    v = -1.0;
	    v = (cj && cju) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = ju;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[ju] -= v;
	  }
	}

	if (jd != j && jd != ju && jd != -1) {
	  if (jd > j) {
	    v = -1.0;
	    v = (cj && cjd) ? -jump : -1.0;
	    /*printf(">> %g\n",v);*/
	    (m->rowind)[ip]  = jd;
	    (m->values.d/*taucs_values*/)[ip]  = v;
	    ip++;
	    (m->values.d/*taucs_values*/)[jp] -= v;
	    D[jd] -= v;
	  }
	}

      }
    }
  }
  taucs_free(D);
  (m->colptr)[N] = ip;

  taucs_printf("taucs_ccs_generate_discontinuous: done, ncols=%d allocated nnz=%d real nnz=%d\n",
	       N,nnz,ip);

  /*taucs_ccs_write_ijv(m,"X.ijv");*/

  return m;
}

#endif /* TAUCS_CORE_DOUBLE */

/*********************************************************/
/*                                                       */
/*********************************************************/
