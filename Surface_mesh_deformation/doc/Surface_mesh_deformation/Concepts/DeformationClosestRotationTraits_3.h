/*!
\ingroup PkgSurfaceMeshDeformationConcepts
\cgalConcept

@brief Concept describing the set of requirements for computing a 3x3 rotation matrix that is close to a given 3x3 matrix, together with basic computations used in the class `CGAL::Surface_mesh_deformation`.
The definition of close depends on the model.
The fact that some basic operations are hidden behind a function is to allow to benefit from optimizations like expression template from libraries used
to implement models of this concept.

\cgalRefines `DefaultConstructible`

\cgalHasModel `CGAL::Deformation_Eigen_closest_rotation_traits_3`
\cgalHasModel `CGAL::Deformation_Eigen_polar_closest_rotation_traits_3`

*/
class DeformationClosestRotationTraits_3{
public:
/// \name Types
/// @{
  /// 3x3 matrix type having a copy constructor and an assignment operator
  typedef unspecified_type Matrix;
  /// 3D vector type having a copy constructor
  typedef unspecified_type Vector;
/// @}

/// \name Operations
/// @{

  /// Equivalent to `result = result + w * (v1*v2^t)`
  void add_scalar_t_vector_t_vector_transpose(Matrix& result, double w, const Vector& v1, const Vector& v2);

  /// Equivalent to `result = result + (w1*m1 + w2*m2) * v`
  void add__scalar_t_matrix_p_scalar_t_matrix__t_vector(Vector& result, double w1, const Matrix& m1, double w2, const Matrix& m2, const Vector& v);

  /// Equivalent to `result = result + w1 * (m1 + m2 + m3) * v`
  void add_scalar_t_matrix_sum_t_vector(Vector& result, double w1, const Matrix& m1, const Matrix& m2, const Matrix& m3, const Vector& v);

  /// Returns the squared norm of `v1 - m*v2`
  double squared_norm_vector_scalar_vector_subs(const Vector& v1, const Matrix& m, const Vector& v2);

  /// Returns an identity matrix
  Matrix identity_matrix();

  /// Returns a matrix initialized with zeros
  Matrix zero_matrix();

  /// Returns the vector (x,y,z)
  Vector vector(double x, double y, double z);

  /// Returns `i`th coefficient of a vector
  double vector_coordinate(const Vector& v, int i);

  /// Computes a rotation matrix close to `m` and places it into `R`
  void compute_close_rotation(const Matrix& m, Matrix& R);

/// @}

};