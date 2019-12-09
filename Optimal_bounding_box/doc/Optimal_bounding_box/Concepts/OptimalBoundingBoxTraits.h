namespace CGAL {

/*!
\ingroup PkgOptimalBoundingBoxConcepts
\cgalConcept

The concept `OptimalBoundingBoxTraits` describes the requirements of the traits class
used in the function `CGAL::optimal_bounding_box()`, and in particular the need for
a 3x3 matrix type with a number of supporting functions (determinant and transposed matrix
computations, etc.).

\cgalRefines `Kernel`

\cgalHasModel `CGAL::Optimal_bounding_box::Optimal_bounding_box_traits`

*/
class OptimalBoundingBoxTraits
{
public:
  /// The field number type; must be a model of the concept `FieldNumberType`.
  typedef unspecified_type                                FT;

  /// A 3x3 matrix type; must be a model of the concept `SvdTraits::Matrix` and support
  /// matrix-matrix and scalar-matrix multiplication, as well as matrix-matrix addition.
  typedef unspecified_type                                Matrix;

  /// Returns the transpose of a matrix.
  Matrix transpose(const Matrix& mat) const;

  /// Returns the determinant of a matrix.
  FT compute_determinant(const Matrix& matrix) const;

  /// Returns the unary matrix Q obtained in the QR decompoisiton of the matrix `A`.
  Matrix get_Q(const Matrix& A) const;
};

} // namespace CGAL
