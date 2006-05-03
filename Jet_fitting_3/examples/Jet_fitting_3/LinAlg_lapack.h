#ifndef _GSL_H_
#define _GSL_H_

#include <stdlib.h>
#include "blaswrap.h"
#include "f2c.h"
extern "C" {
#include "clapack.h"
}
//#include <math.h>


////////////////////////class gsl_Matrix/////////////////////
//in lapack, matrices are one-dimensional arrays 
// and elements are column-major ordered
class gsl_Matrix{
protected:
  double* m_matrix;
public:
  size_t nb_rows;
  size_t nb_columns;
  //contructor
  // initializes all the elements of the matrix to zero.
  gsl_Matrix(size_t n1, size_t n2) { 
    m_matrix = (double*) calloc (n1*n2, sizeof(double)); 
    nb_rows = n1;
    nb_columns = n2;
  }

  //access
  const double* matrix() const { return m_matrix; }
  double* matrix() { return m_matrix; }

  void set_elt(size_t i, size_t j, const double value) { m_matrix[j*nb_rows+i] = value; }
  double get_elt(size_t i, size_t j) { return m_matrix[j*nb_rows+i]; }
}; 

////////////////////////class GSL/////////////////////
class GSL{
public:
  //  typedef double FT;
  typedef gsl_Matrix Matrix;
  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order, 
  // eigen vectors are sorted in accordance.
  static
    void eigen_symm_algo(Matrix& S, double* eval, Matrix& evec);
  //solve MX=B using SVD and give the condition number of M
  //M replaced by U after gsl_linalg_SV_decomp()
  //The diagonal and lower triangular part of M are destroyed during the
  //computation, but the strict upper triangular part is not referenced.
  static
    void solve_ls_svd_algo(Matrix& M, double* B, double &cond_nb);
};

void GSL::eigen_symm_algo(Matrix& S, double* eval, Matrix& evec)
{
  char jobz = 'V';
  char range = 'A';
  char uplo = 'U';
  integer n = S.nb_rows,
    lda = S.nb_rows,
    il = 0, iu = 0,
    m = n, 
    lwork = 26*n,
    liwork = 10*n,
    info;
  double vl=0, vu=0, abstol = 0;
/*   //c style */
/*   integer* isuppz = (integer*) malloc(n*sizeof(integer)); */
/*   double* work = (double*) malloc(lwork*sizeof(double)); */
/*   integer* iwork = (integer*) malloc(liwork*sizeof(integer)); */
  //c++ style
   integer* isuppz;
  isuppz = new integer [n];
  double* work = new double [lwork];
  integer* iwork = new integer [liwork]; 

  dsyevr_(&jobz, &range, &uplo, &n, S.matrix(), &lda, &vl, &vu, &il, &iu, &abstol, &m,
	  eval, evec.matrix(), &n, isuppz, work, &lwork, iwork, &liwork, &info);
   
/*   //clean up  */
/*   free(isuppz); */
/*   free(work); */
/*   free(iwork); */
/*    //C++ style */
   delete [] isuppz;
   delete [] work;
   delete [] iwork;
}

void GSL::solve_ls_svd_algo(Matrix& M, double* B, double &cond_nb)
{
  integer m = M.nb_rows,
    n = M.nb_columns,
    nrhs = 1,
    lda = m,
    ldb = m,
    rank,
    lwork = 5*m,
    info;
/*   //c style */
/*   double* sing_values = (double*) malloc(n*sizeof(double)); */
/*   double* work = (double*) malloc(lwork*sizeof(double)); */
  //c++ style
  double* sing_values = new double [n];
  double* work = new double [lwork];
 
  double rcond = -1;

  dgelss_(&m, &n, &nrhs, M.matrix(), &lda, B, &ldb, sing_values, 
	  &rcond, &rank, work, &lwork, &info);
  assert(info==0);

  cond_nb = sing_values[0]/sing_values[n-1];
  
/*   //clean up  */
/*   free(sing_values); */
/*   free(work); */
//c++ style
   delete [] sing_values;
  delete [] work;
}

#endif
