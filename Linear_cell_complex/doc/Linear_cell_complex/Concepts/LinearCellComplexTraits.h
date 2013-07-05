
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

Required types and functors for the `LinearCellComplexTraits` concept. This
geometric traits concept is used in the \ref CGAL::Linear_cell_complex "Linear_cell_complex"
class.

\cgalHasModel \ref CGAL::Linear_cell_complex_traits "CGAL::Linear_cell_complex_traits<d,K>"

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`

*/

class LinearCellComplexTraits {
public:

/// \name Constants
/// @{

/*! The ambient dimension (must be > 1).
*/
static unsigned int ambient_dimension;

/// @}

/// \name Types
/// @{

/*!
a number type that is a model of FieldNumberType.
*/
typedef unspecified_type FT;

/*!
point type.
*/
typedef unspecified_type Point;

/*!
vector type.
*/
typedef unspecified_type Vector;

/// @}

/// \name Constructions
/// @{

/*!
Functor that provides \ref LinearCellComplexTraits::Point "Point " `operator() (const` \ref Point " Point"`& p, const` \ref LinearCellComplexTraits::Vector " Vector"` & v)`,
which constructs the translation of point `p` by vector `v`, and
\ref LinearCellComplexTraits::Point "Point " `operator() (const CGAL::Origin&, const ` \ref LinearCellComplexTraits::Vector " Vector"& v)`,
which constructs the translation of a point at the origin by vector `v`
(used in \ref CGAL::Linear_cell_complex::barycenter "Linear_cell_complex::barycenter").
*/
typedef unspecified_type Construct_translated_point;

/*!
Functor that provides \ref LinearCellComplexTraits::Vector " Vector " `operator() (const ` \ref Point "Point"`& p1, const ` \ref Point " Point"`& p2)`
which constructs a vector as the difference of points `p2-p1`, and
\ref LinearCellComplexTraits::Vector " Vector " `operator() (const CGAL::Origin&, const ` \ref Point " Point"`& p)`
which constructs a vector as the difference of point `p` and a point at the origin
(used in \ref CGAL::Linear_cell_complex::barycenter "Linear_cell_complex::barycenter"
and `CGAL::import_from_plane_graph`).
*/
typedef unspecified_type Construct_vector;

/*!
Functor that provides \ref LinearCellComplexTraits::Vector " Vector " `operator() (const` \ref LinearCellComplexTraits::Vector " Vector"`& v1, const` \ref LinearCellComplexTraits::Vector " Vector"`& v2)`
which constructs a vector as the sum of vectors `v1+v2`
(used in \ref CGAL::Linear_cell_complex::barycenter "Linear_cell_complex::barycenter",
`CGAL::compute_normal_of_cell_0`
and `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_sum_of_vectors;

/*!
Functor that provides \ref LinearCellComplexTraits::Vector " Vector " `operator() (const` \ref LinearCellComplexTraits::Vector " Vector"`& v, ` \ref LinearCellComplexTraits::FT "FT" `scale)`
which constructs a vector equal to vector `v` scaled by `scale` factor
(used in \ref CGAL::Linear_cell_complex::barycenter "Linear_cell_complex::barycenter",
`CGAL::compute_normal_of_cell_0` and `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_scaled_vector;

/*!
Functor that provides \ref LinearCellComplexTraits::Point "Point " `operator() (const ` \ref Point "Point"`& p1, const ` \ref Point "Point"`& p2)`
which constructs the midpoint of points `p1` and `p2`
(used in \ref CGAL::Linear_cell_complex::barycenter "Linear_cell_complex::barycenter").
*/
typedef unspecified_type Construct_midpoint;

/// @}

/// \name
/// If `ambient_dimension==2`
/// @{

/*!
a model of \ref Kernel::Direction_2 "Direction_2".
*/
typedef unspecified_type Direction_2;

/*!
a model of \ref Kernel::ConstructDirection_2 "ConstructDirection_2" (used in `CGAL::import_from_plane_graph`).
*/
typedef unspecified_type Construct_direction_2;

/// @}

/// \name
/// If `ambient_dimension==3`
/// @{

/*!
a model of \ref Kernel::ConstructNormal_3 "ConstructNormal_3" (used in `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Construct_normal_3;

/*!
a model of \ref Kernel::Collinear_3 "Collinear_3" (used in `CGAL::compute_normal_of_cell_2`).
*/
typedef unspecified_type Collinear_3;

/// @}

}; /* end LinearCellComplexTraits */

