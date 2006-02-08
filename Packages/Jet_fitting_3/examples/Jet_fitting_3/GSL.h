#ifndef _GSL_H_
#define _GSL_H_

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

////////////////////////class gsl_Vector/////////////////////
class gsl_Vector{
protected:
  gsl_vector* m_vector;
public:
  //contructor 
  // initializes all the elements of the vector to zero
  gsl_Vector(size_t n) { m_vector = gsl_vector_calloc(n); }
  //access
  double const operator[](size_t i) const { return gsl_vector_get(m_vector,
								  i); }
  //use as left value v[i]=10 
  double& operator[](size_t i) { return *gsl_vector_ptr(m_vector, i);}
  const gsl_vector* vector() const { return m_vector; }
  gsl_vector* vector() { return m_vector; }
};

////////////////////////class gsl_Matrix/////////////////////
class gsl_Matrix{
protected:
  gsl_matrix* m_matrix;
public:
  //contructor
  // initializes all the elements of the matrix to zero.
  gsl_Matrix(size_t n1, size_t n2) { m_matrix = gsl_matrix_calloc (n1, n2); }

  //class Row, to define the usual double operator [][] for matrices
 class Row{
    gsl_matrix* matrix;
    size_t row;
  public:
    Row(gsl_matrix* _matrix, size_t _row) : 
      matrix(_matrix), row(_row) {}
    double const operator[](size_t column) const 
      { return gsl_matrix_get(this->matrix, row, column);}
    double& operator[](size_t column)
      { return *gsl_matrix_ptr(this->matrix, row, column);}
  };//END class Row

  Row operator[](size_t _row) { return Row(m_matrix, _row);}

  //access
  const gsl_matrix* matrix() const { return m_matrix; }
  gsl_matrix* matrix() { return m_matrix; }
}; 

////////////////////////class GSL/////////////////////
class GSL{
public:
  typedef gsl_Vector Vector;
  typedef gsl_Matrix Matrix;
  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  static
  void eigen_symm_algo(Matrix& S, Vector& eval, Matrix& evec);
  //solve MX=B using SVD and give the condition number of M
  //M replaced by U after gsl_linalg_SV_decomp()
  //The diagonal and lower triangular part of M are destroyed during the
  //computation, but the strict upper triangular part is not referenced.
  static
  void solve_ls_svd_algo(Matrix& M, Vector& X, const Vector& B,
			 double &cond_nb);
};

void GSL::eigen_symm_algo(Matrix& S, Vector& eval, Matrix& evec)
{
  const size_t common_size = S.matrix()->size1;
  assert( S.matrix()->size2 == common_size
	  && eval.vector()->size == common_size
	  && evec.matrix()->size1 == common_size
	  && evec.matrix()->size2 == common_size
	  );
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (common_size);
  gsl_eigen_symmv (S.matrix(), eval.vector(), evec.matrix(), w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (eval.vector(), evec.matrix(),
			GSL_EIGEN_SORT_VAL_DESC); 
}

void GSL::solve_ls_svd_algo(Matrix& M, Vector& X, const Vector& B,
			    double &cond_nb)
{
  const size_t nb_lines = M.matrix()->size1;
  assert( B.vector()->size == nb_lines );
  const size_t nb_columns = M.matrix()->size2;
  assert( X.vector()->size == nb_columns );

  gsl_matrix * V = gsl_matrix_alloc(nb_columns,nb_columns);
  gsl_vector * S = gsl_vector_alloc(nb_columns);
  gsl_vector * work = gsl_vector_alloc(nb_columns);

  gsl_linalg_SV_decomp (M.matrix(), V, S, work);
  gsl_linalg_SV_solve (M.matrix(), V, S, B.vector(), X.vector());

  cond_nb = gsl_vector_get(S,0)/gsl_vector_get(S,nb_columns-1);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
}

#endif
