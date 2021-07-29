namespace CGAL {
namespace Barycentric_coordinates {

/*!
\ingroup PkgBarycentricCoordinates3RefConcepts
\cgalConcept

A concept that describes the set of requirements of the template parameter
`GeomTraits` used to parameterize all classes and functions with 3D barycentric
coordinates from the namespace `CGAL::Barycentric_coordinates`.

\cgalHasModel
- All models of `Kernel`
*/
class BarycentricTraits_3 {

public:

/// \name Types
/// @{

/*!
  A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/// @}

/// \name 3D Geometric Objects
/// @{

/*!
  A model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
  A model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
  A model of `Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/// @}

/// \name 3D Generalized Constructions
/// @{

/*!
  A construction object that must provide the function operator:

  `FT operator(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns the signed area of the triangle defined by the points `p`, `q`, and `r`.
*/
typedef unspecified_type Compute_area_2;

/*!
  A construction object that must provide the function operator:

  `FT operator(const Point_3& p0, const Point_3& p1, const Point_3& p2, const Point_3& p3)`
  
  that returns the signed volume of the tetrahedron defined by the four points `p0`, `p1`, `p2`, and `p3`.
*/
typedef unspecified_type Compute_volume_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(const Vector_3& u, const Vector_3& v)`
  
  that returns an approximation of the angle between `u` and `v`.
  The angle is given in degrees.
*/
typedef unspecified_type Compute_approximate_angle_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(const Vector_3& v, const Vector_3& w)`
  
  that returns the scalar (inner) product of the two vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(const Vector_3& u, const Vector_3& v, const Vector_3& w)`
  
  that returns the determinant of the three vectors `u`, `v` and `w`.
*/
typedef unspecified_type Compute_determinant_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(const FT& value)`
  
  that returns the square root of value.
*/
typedef unspecified_type Sqrt;

/*!
  A construction object that must provide the function operator:

  `Vector_3 operator(const Point_3& p, const Point_3& q)`
  
  that returns the vector `q` - `p`.
*/
typedef unspecified_type Construct_vector_3;

/*!
  A construction object that must provide the function operator:

  `Vector_3 operator(const Vector_3& u, const Vector_3& v)`
  
  that returns the cross product between `u` and `v`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/// @}

/// \name 2D Generalized Predicates
/// @{


/// @}

};

} // namespace Barycentric_coordinates
} // namespace CGAL
