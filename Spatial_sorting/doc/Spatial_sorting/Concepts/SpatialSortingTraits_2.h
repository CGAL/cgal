
/*!
\ingroup PkgSpatialSortingConcepts
\cgalConcept

All 2D spatial sorting algorithms provided in \cgal are parameterized
by a traits class `Traits`, which defines the
primitives (objects and predicates) that the sorting algorithms use.
`SpatialSortingTraits_2` defines the complete set of primitives required in these
functions and functors.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModels{CGAL::Spatial_sort_traits_adapter_2<Base_traits,PointPropertyMap>}
\cgalHasModelsEnd

*/

class SpatialSortingTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type on which the sorting algorithms operate.
*/
typedef unspecified_type Point_2;

/*!
Functor object type returning the \f$ x\f$ coordinate of a `Point_2`.
Must provide
`FT operator()(Point_2 p)` where `FT` can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_x_2;

/*!
Functor object type returning the \f$ y\f$ coordinate of a `Point_2`.
Must provide
`FT operator()(Point_2 p)` where `FT` can be used as
argument of `CGAL::to_double`.

*/
typedef unspecified_type Compute_y_2;

/*!
Binary predicate object type comparing `Point_2`s
along the \f$ x\f$ coordinate.
Must provide
`bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p_x < q_x\f$,
where \f$ p_x\f$ and \f$ q_x\f$ denote \f$ x\f$ coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_x_2;

/*!
Binary predicate object type comparing `Point_2`s
along the \f$ y\f$ coordinate.
Must provide
`bool operator()(Point_2 p, Point_2 q)` where `true`
is returned iff \f$ p_y < q_y\f$,
where \f$ p_y\f$ and \f$ q_y\f$ denote \f$ y\f$ coordinate of point \f$ p\f$ and \f$ q\f$,
respectively.

*/
typedef unspecified_type Less_y_2;

/// @}

/// \name Creation
/// Only a copy constructor is required.
/// @{

/*!

*/
SpatialSortingTraits_2(const SpatialSortingTraits_2& t);

/// @}

/// \name Operations
/// The following member functions to create instances of the above predicate object types must exist.
/// @{

/*!

*/
Compute_x_2 compute_x_2_object();

/*!

*/
Compute_y_2 compute_y_2_object();

/*!

*/
Less_x_2 less_x_2_object();

/*!

*/
Less_y_2 less_y_2_object();

/// @}

}; /* end SpatialSortingTraits_2 */

