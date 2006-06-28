#ifndef CGAL_LAPACK_H
#define CGAL_LAPACK_H

#include <stdlib.h>
#include "blaswrap.h"
#include "f2c.h"
extern "C" {
#include "clapack.h"
}

namespace CGAL {
////////////////////////class Lapack_matrix/////////////////////
//in clapack, matrices are one-dimensional arrays and elements are
//column-major ordered. This class is a wrapper defining set and get
//in the usual way with line and column indices.
class Lapack_matrix{
protected:
  double* m_matrix;
public:
  size_t nb_rows;
  size_t nb_columns;
  //contructor
  // initializes all the elements of the matrix to zero.
  Lapack_matrix(size_t n1, size_t n2) { 
    m_matrix = (double*) calloc (n1*n2, sizeof(double)); 
    nb_rows = n1;
    nb_columns = n2;
  }

  //destructor
  //  ~Lapack_matrix();// { free m_matrix; }
  
  //access
  const double* matrix() const { return m_matrix; }
  double* matrix() { return m_matrix; }

  void set_elt(size_t i, size_t j, const double value) { m_matrix[j*nb_rows+i] = value; }
  double get_elt(size_t i, size_t j) { return m_matrix[j*nb_rows+i]; }
}; 

//Lapack_matrix::~Lapack_matrix() {  free m_matrix; } bug!

////////////////////////class Lapack/////////////////////
class Lapack{
public:
  typedef Lapack_matrix Matrix;
 //solve MX=B using SVD and give the condition number of M
  static
    void solve_ls_svd_algo(Matrix& M, double* B, double &cond_nb);
};

void Lapack::solve_ls_svd_algo(Matrix& M, double* B, double &cond_nb)
{
  integer m = M.nb_rows,
    n = M.nb_columns,
    nrhs = 1,
    lda = m,
    ldb = m,
    rank,
    lwork = 5*m,
    info;
  //c style
  double* sing_values = (double*) malloc(n*sizeof(double));
  double* work = (double*) malloc(lwork*sizeof(double));

  double rcond = -1;

  dgelss_(&m, &n, &nrhs, M.matrix(), &lda, B, &ldb, sing_values, 
	  &rcond, &rank, work, &lwork, &info);
  assert(info==0);

  cond_nb = sing_values[0]/sing_values[n-1];
  
  //clean up 
  free(sing_values);
  free(work);
}

} // namespace CGAL

#endif // CGAL_LAPACK_H
