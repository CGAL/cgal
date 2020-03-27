/*!
\ingroup PkgOptimalBoundingBoxConcepts
\cgalConcept

The concept `OrientedBoundingBoxTraits_3` describes the requirements of the traits class
used in the function `CGAL::oriented_bounding_box()`, and in particular the need for
a 3x3 matrix type.

\cgalRefines `Kernel`

\cgalHasModel `CGAL::Oriented_bounding_box_traits_3`

*/
class OrientedBoundingBoxTraits_3
{
public:
  /// The field number type; must be a model of the concept `FieldNumberType`
  typedef unspecified_type                                FT;

  /// The 3D affine transformation type; the template parameter `K` must be a model of `Kernel`
  /// and be compatible with the type `Point_3`
  typedef CGAL::Aff_transformation_3<K>                   Aff_transformation_3;

  /// A 3x3 matrix type; model of the concept `SvdTraits::Matrix` and which supports
  /// matrix-matrix and scalar-matrix multiplication, as well as matrix-matrix addition
  typedef unspecified_type                                Matrix;

  /// A 3 dimensional vector type; model of the concept `SvdTraits::Vector` and which supports
  /// matrix-vector and multiplication
  typedef unspecified_type                                Vector;

  /// Returns the unitary matrix `Q` obtained in the QR-decomposition of the matrix `m`
  Matrix get_Q(const Matrix& m) const;
};
