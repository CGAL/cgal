/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

@brief Concept describing the set of requirements for a SVD factorization of a 3x3 matrix together with basic computations used in the class `CGAL::Deform_mesh`.
The fact that some basic operations are hidden behind a function is to allow to benefit from optimizations like expression template from libraries used
to implement a model of this concept.

\cgalHasModel `CGAL::Deformation_Eigen_Svd_traits_3`

*/
class DeformationSvdTraits_3{
public:
/// \name Types 
/// @{
  /// <I>3x3</I> matrix type with copy constructor and assignment operator
  typedef Hidden_type Matrix;
  /// <I>3x1</I> vector type with copy constructor
  typedef Hidden_type Vector;
  /// SVD factorization solver for DeformationSvdTraits_3::Matrix
  typedef Hidden_type Solver;
/// @} 

/// \name Creation
/// @{
/// Default constructor.
  DeformationSvdTraits_3();
/// @}

/// \name Operations 
/// @{

  /// Equivalent to `result = m1*m2^t`
  void matrix_matrix_transpose_mult(Matrix& result, const Matrix& m1, const Matrix& m2);
  
  /// Equivalent to `result += w * (v1*v2^t)`
  void scalar_vector_vector_transpose_mult(Matrix& result, double w, const Vector& v1, const Vector& v2);
  
  /// Equivalent to `result += (w1*m1 + w2*m2) * v`
  void scalar_matrix_scalar_matrix_vector_mult(Vector& result, double w1, const Matrix& m1, double w2, const Matrix& m2, const Vector& v);
  
  /// Equivalent to `result += w1 * (m1 + m2 + m3) * v`
  void scalar_mult_with_matrix_sum(Vector& result, double w1, const Matrix& m1, const Matrix& m2, const Matrix& m3, const Vector& v);
  
  /// Returns the squared norm of `v1 - m*v2`
  double squared_norm_vector_scalar_vector_subs(const Vector& v1, const Matrix& m, const Vector& v2);

  /// Negates a column of a matrix
  void negate_column(Matrix& result, int i);

  /// Returns an identity matrix
  Matrix identity_matrix();
  
  /// Returns a zero initialized matrix
  Matrix zero_matrix();

  /// Returns a vector initialized with parameters
  Vector vector(double x, double y, double z);
  
  /// Returns `i`th coefficient of a vector
  double vector_coeff(const Vector& v, int i);
  
  /// Returns determinant of `m`
  double determinant(const Matrix& m);
    
  /// Computes the singular value decomposition and returns the solver
  Solver compute_svd(const Matrix& m);
  
  /// Returns the diagonal index of smallest singular value 	
  int get_smallest_singular_value_index(const Solver& solver);
  
  /// Gets matrix U from solver
  const Matrix& get_matrixU(const Solver& solver);

  /// Gets matrix V from solver
  const Matrix& get_matrixV(const Solver& solver);
/// @}
 
};