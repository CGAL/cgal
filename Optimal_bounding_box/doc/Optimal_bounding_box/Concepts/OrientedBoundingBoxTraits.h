namespace CGAL {

/*!
\ingroup PkgOrientedBoundingBoxConcepts
\cgalConcept

The concept `OrientedBoundingBoxTraits` describes the requirements of the traits class
used in the function `CGAL::oriented_bounding_box()`, and in particular the need for
a 3x3 matrix type with a number of supporting functions (determinant and transposed matrix
computations, etc.).

\cgalHasModel `CGAL::Oriented_bounding_box_traits`

*/
class OrientedBoundingBoxTraits
{
public:
  /// The field number type; must be a model of the concept `FieldNumberType`.
  typedef unspecified_type                                FT;

  /// The 3D point type; must be model of `Point_3`
  typedef unspecified_type                                Point_3;

  /// The 3D affine transformation type; the template parameter `K` must be a model of `Kernel`
  /// and be compatible with the type `Point_3`.
  typedef CGAL::Aff_transformation_3<K>                   Aff_transformation_3;

  /// A construction object that must provide the function operator:
  /// `CGAL::Bbox_3 operator()(const Point_3&)`,
  /// which returns an axis-aligned bounding that contains the point
  typedef unspecified_type                               Construct_bbox_3;

  /// A 3x3 matrix type; model of the concept `SvdTraits::Matrix` and which supports
  /// matrix-matrix and scalar-matrix multiplication, as well as matrix-matrix addition.
  typedef unspecified_type                                Matrix;

  /// Returns the transpose of the matrix `m`.
  Matrix transpose(const Matrix& m) const;

  /// Returns the determinant of the matrix `m`.
  FT compute_determinant(const Matrix& m) const;

  /// Returns the unary matrix `Q` obtained in the QR-decomposition of the matrix `m`.
  Matrix get_Q(const Matrix& m) const;
};

} // namespace CGAL
