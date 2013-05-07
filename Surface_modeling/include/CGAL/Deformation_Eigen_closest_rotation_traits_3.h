#ifndef CGAL_DEFORMATION_EIGEN_CLOSEST_ROTATION_TRAITS_3_H
#define CGAL_DEFORMATION_EIGEN_CLOSEST_ROTATION_TRAITS_3_H

#include <Eigen/Eigen>
#include <Eigen/SVD>

namespace CGAL {
/// \ingroup PkgSurfaceModeling
/// A class to compute the closest rotation in Frobenius norm to a 3x3 Matrix using `Eigen` library.
/// The internal computation relies on `Eigen::JacobiSVD<>` solver.
///
/// \cgalModels `DeformationClosestRotationTraits_3`
class Deformation_Eigen_closest_rotation_traits_3{
public:

  /// \cond SKIP_FROM_MANUAL
  typedef Eigen::Matrix3d Matrix;
  typedef Eigen::Vector3d Vector;

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
    
  /// Computes the closest rotation to `m` and places it into `R`
  void compute_close_rotation(const Matrix& m, Matrix& R)
  {
    Eigen::JacobiSVD<Eigen::Matrix3d> solver;
    solver.compute( m, Eigen::ComputeFullU | Eigen::ComputeFullV );

    const Matrix& u = solver.matrixU(); const Matrix& v = solver.matrixV();
    R = v * u.transpose();

    if( R.determinant() < 0 ) {
      Matrix u_copy = u;
      u_copy.col(2) *= -1;        // singular values sorted ascendingly  
      R = v * u_copy.transpose(); // re-extract rotation matrix
    }
  }

  /// \endcond

};

}//namespace CGAL
#endif // CGAL_DEFORMATION_EIGEN_CLOSEST_ROTATION_TRAITS_3_H