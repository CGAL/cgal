#ifndef CGAL_EIGEN_SVD_HELPER_3_H
#define CGAL_EIGEN_SVD_HELPER_3_H

#include <Eigen/Eigen>
#include <Eigen/SVD>

namespace CGAL {
/// \ingroup PkgSurfaceModeling
/// Class providing functionality to compute SVD factorization of a 3x3 Matrix . `Eigen::JacobiSVD<>` is used as internal SVD solver.
///
/// \cgalModels `SVDHelper_3`
class Eigen_SVD_helper_3{
public:

  typedef Eigen::Matrix3d Matrix;
  typedef Eigen::Vector3d Vector;
  typedef Eigen::JacobiSVD<Eigen::Matrix3d> Solver;
  
  /// Equivalent to `result = m1*m2^t`
  void matrix_matrix_transpose_mult(Matrix& result, const Matrix& m1, const Matrix& m2) 
  { 
    result = m1*m2.transpose();
  }

  /// Equivalent to `result += w * (v1*v2^t)`
  void scalar_vector_vector_transpose_mult(Matrix& result, double w, const Vector& v1, const Vector& v2) 
  {
    result += w * (v1*v2.transpose());
  }
  
  /// Equivalent to `result += (w1*m1 + w2*m2) * v`
  void scalar_matrix_scalar_matrix_vector_mult(Vector& result, double w1, const Matrix& m1, double w2, const Matrix& m2, const Vector& v) 
  {
    result += (w1*m1 + w2*m2) * v;
  }
  
  /// Equivalent to `result += w * (m1 + m2 + m3) * v`
  void scalar_mult_with_matrix_sum(Vector& result, double w, const Matrix& m1, const Matrix& m2, const Matrix& m3, const Vector& v) 
  { 
    result += w * (m1 + m2 + m3) * v;
  }
  
  /// Returns the squared norm of `v1 - m*v2`
  double squared_norm_vector_scalar_vector_subs(const Vector& v1, const Matrix& m, const Vector& v2) 
  { 
    return (v1 - m*v2).squaredNorm(); 
  } 

  /// Negates column `i` of matrix `result`
  void negate_column(Matrix& result, int i) 
  {
    result.col(i) *= -1;  
  }

  /// Returns an identity matrix
  Matrix identity_matrix() 
  {
    return Matrix().setIdentity();
  }
  
  /// Returns a zero initialized matrix
  Matrix zero_matrix() 
  {
    return Matrix().setZero();
  }

  /// Returns vector initialized with parameters
  Vector vector(double x, double y, double z)
  {
    return Vector(x, y, z);
  }
  
  /// Returns a coefficient of a vector
  double vector_coeff(const Vector& v, int i)
  {
    return v(i);
  }
  
  /// Returns determinant of m
  double determinant(const Matrix& m) 
  {
    return m.determinant();
  }
    
  /// Computes the singular value decomposition and returns the solver
  Solver compute_svd(const Matrix& m)
  {
    return Solver().compute( m, Eigen::ComputeFullU | Eigen::ComputeFullV );
  }
	
  /// Returns the diagonal index of smallest singular value 	
  int get_smallest_singular_value_index(const Solver& /*solver*/)
  {
    return 2; // singular values are always sorted in decreasing order so use column 2
  }

  /// Gets matrix U from solver
  const Matrix& get_matrixU(const Solver& solver) 
  {
    return solver.matrixU();
  }

  /// Gets matrix V from solver
  const Matrix& get_matrixV(const Solver& solver) 
  {
    return solver.matrixV();
  }
};

}//namespace CGAL
#endif // CGAL_EIGEN_SVD_HELPER_3_H