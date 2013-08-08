/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

@brief Concept describing the set of requirements for computing a close rotation to a 3x3 matrix together with basic computations used in the class `CGAL::Deform_mesh`.
The fact that some basic operations are hidden behind a function is to allow to benefit from optimizations like expression template from libraries used
to implement models of this concept.

\cgalHasModel `CGAL::Deformation_Eigen_closest_rotation_traits_3`
\cgalHasModel `CGAL::Deformation_Eigen_polar_closest_rotation_traits_3`

*/
class DeformationClosestRotationTraits_3{
public:
/// \name Types 
/// @{
  /// <I>3x3</I> matrix type having a copy constructor and an assignment operator
  typedef Hidden_type Matrix;
  /// <I>3D</I> vector type having a copy constructor
  typedef Hidden_type Vector;
/// @} 

/// \name Creation
/// @{
/// Default constructor.
  DeformationClosestRotationTraits_3();
/// @}

/// \name Operations 
/// @{
  
  /// Equivalent to `result += w * (v1*v2^t)`
  void scalar_vector_vector_transpose_mult(Matrix& result, double w, const Vector& v1, const Vector& v2);
  
  /// Equivalent to `result += (w1*m1 + w2*m2) * v`
  void scalar_matrix_scalar_matrix_vector_mult(Vector& result, double w1, const Matrix& m1, double w2, const Matrix& m2, const Vector& v);
  
  /// Equivalent to `result += w1 * (m1 + m2 + m3) * v`
  void scalar_mult_with_matrix_sum(Vector& result, double w1, const Matrix& m1, const Matrix& m2, const Matrix& m3, const Vector& v);
  
  /// Returns the squared norm of `v1 - m*v2`
  double squared_norm_vector_scalar_vector_subs(const Vector& v1, const Matrix& m, const Vector& v2);

  /// Returns an identity matrix
  Matrix identity_matrix();
  
  /// Returns a matrix initialized with zero
  Matrix zero_matrix();

  /// Returns the vector (x,y,z)
  Vector vector(double x, double y, double z);
  
  /// Returns `i`th coefficient of a vector
  double vector_coeff(const Vector& v, int i);
    
  /// Computes a rotation matrix close to `m` and places it into `R`
  void compute_close_rotation(const Matrix& m, Matrix& R);
  
/// @}
 
};