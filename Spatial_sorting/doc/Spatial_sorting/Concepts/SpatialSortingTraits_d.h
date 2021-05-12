
/*!
\ingroup PkgSpatialSortingConcepts
\cgalConcept

All higher dimensional spatial sorting algorithms provided in \cgal are parameterized
by a traits class `Traits`, which defines the
primitives (objects and predicates) that the sorting algorithms use.
`SpatialSortingTraits_d` defines the complete set of primitives required in these
functions and functors.

\cgalHasModel Any \cgal `d`-dimensional kernel.
\cgalHasModel `CGAL::Spatial_sort_traits_adapter_d<Base_traits,PointPropertyMap>`

*/

class SpatialSortingTraits_d {
public:

/// \name Types
/// @{

/*!
The point type on which the sorting algorithms operate.

*/
typedef unspecified_type Point_d;

/*!
Functor object type returning the dimension of a `Point_d`.
Must provide
`int operator()(Point_d p)` returning the dimension of \f$ p\f$.

*/
typedef unspecified_type Point_dimension_d;

/*!
Functor object type returning the coordinates of a `Point_d`.
Must provide
`FT operator()(Point_d p, int i)` returning the \f$ i\f$th
coordinate of \f$ p\f$. `FT` is a type that can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_coordinate_d;

/*!
Binary predicate object type comparing `Point_d`s
along some coordinate.
Must provide
`bool operator()(Point_d p, Point_d q, int i)` where `true`
is returned iff \f$ p_i < q_i\f$,
where \f$ p_i\f$ and \f$ q_i\f$ denote \f$ i\f$th coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_coordinate_d;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
SpatialSortingTraits_d(const SpatialSortingTraits_d& t);

/// @}

/// \name Operations
/// The following member functions to create instances of the above predicate object types must exist.
/// @{

/*!

*/
Point_dimension_d point_dimension_d_object();

/*!

*/
Compute_coordinate_d compute_coordinate_d_object();

/*!

*/
Less_coordinate_d less_coordinate_d_object();

/// @}

}; /* end SpatialSortingTraits_d */

