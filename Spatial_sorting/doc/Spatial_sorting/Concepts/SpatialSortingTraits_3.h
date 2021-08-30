
/*!
\ingroup PkgSpatialSortingConcepts
\cgalConcept

All 3D spatial sorting algorithms, including sorting on the sphere, provided in \cgal are parameterized
by a traits class `Traits`, which defines the
primitives (objects and predicates) that the sorting algorithms use.
`SpatialSortingTraits_3` defines the complete set of primitives required in these
functions and functors.

\cgalHasModel Any \cgal kernel.
\cgalHasModel `CGAL::Spatial_sort_traits_adapter_3<Base_traits,PointPropertyMap>`

*/

class SpatialSortingTraits_3 {
public:

/// \name Types
/// @{

/*!
The point type on which the sorting algorithms operate.
*/
typedef unspecified_type Point_3;

/*!
Functor object type returning the \f$ x\f$ coordinate of a `Point_3`.
Must provide
`FT operator()(Point_3 p)` where `FT` can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_x_3;

/*!
Functor object type returning the \f$ y\f$ coordinate of a `Point_3`.
Must provide
`FT operator()(Point_3 p)` where `FT` can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_y_3;

/*!
Functor object type returning the \f$ z\f$ coordinate of a `Point_3`.
Must provide
`FT operator()(Point_3 p)` where `FT` can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_z_3;

/*!
Binary predicate object type comparing `Point_3`s
along the \f$ x\f$ coordinate.
Must provide
`bool operator()(Point_3 p, Point_3 q)` where `true`
is returned iff \f$ p_x < q_x\f$,
where \f$ p_x\f$ and \f$ q_x\f$ denote \f$ x\f$ coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_x_3;

/*!
Binary predicate object type comparing `Point_3`s
along the \f$ y\f$ coordinate.
Must provide
`bool operator()(Point_3 p, Point_3 q)` where `true`
is returned iff \f$ p_y < q_y\f$,
where \f$ p_y\f$ and \f$ q_y\f$ denote \f$ y\f$ coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_y_3;

/*!
Binary predicate object type comparing `Point_3`s
along the \f$ z\f$ coordinate.
Must provide
`bool operator()(Point_3 p, Point_3 q)` where `true`
is returned iff \f$ p_z < q_z\f$,
where \f$ p_z\f$ and \f$ q_z\f$ denote \f$ z\f$ coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_z_3;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
SpatialSortingTraits_3(const SpatialSortingTraits_3& t);

/// @}

/// \name Operations
/// The following member functions to create instances of the above predicate object types must exist.
/// @{

/*!

*/
Compute_x_3 compute_x_3_object();

/*!

*/
Compute_y_3 compute_y_3_object();

/*!

*/
Compute_z_3 compute_z_3_object();

/*!

*/
Less_x_3 less_x_3_object();

/*!

*/
Less_y_3 less_y_3_object();

/*!

*/
Less_z_3 less_z_3_object();

/// @}

}; /* end SpatialSortingTraits_3 */

