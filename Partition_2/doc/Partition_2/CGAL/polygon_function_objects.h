namespace CGAL {

/*!
\ingroup PkgPartition2FunctionObjectClasses

Function object class for testing if a sequence of points represents
a convex polygon or not.

\cgalModels `PolygonIsValid`

\sa `CGAL::convex_partition_is_valid_2()`
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`

\cgalHeading{Implementation}

This test requires \f$ O(n)\f$ time for a polygon with \f$ n\f$ vertices.

*/
template< typename Traits >
class Is_convex_2 {
public:

/// \name Creation
/// @{

/*!
`Traits` satisfies the
requirements of the function `is_convex_2()`
*/
Is_convex_2(const Traits& t);

/// @}

/// \name Operations
/// @{

/*!

returns `true` iff the points of type `Traits::Point_2`
in the range [`first`,`beyond`) define a convex polygon.

*/
template<class InputIterator>
bool operator()(InputIterator first, InputIterator beyond);

/// @}

}; /* end Is_convex_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPartition2FunctionObjectClasses

Function object class that indicates all sequences of points are valid.

\cgalModels `PolygonIsValid`

\sa `CGAL::partition_is_valid_2()`
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`

\cgalHeading{Implementation}

This test requires \f$ O(1)\f$ time.

*/
template< typename Traits >
class Is_vacuously_valid {
public:

/// \name Creation
/// @{

/*!

*/
Is_vacuously_valid(const Traits& t);

/// @}

/// \name Operations
/// @{

/*!

returns `true`.

*/
template<class InputIterator>
bool operator()(InputIterator first, InputIterator beyond);

/// @}

}; /* end Is_vacuously_valid */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgPartition2FunctionObjectClasses

Function object class that tests whether a sequence of points represents
a \f$ y\f$-monotone polygon or not.

\cgalModels `PolygonIsValid`

\sa `CGAL::convex_partition_is_valid_2()`
\sa `CGAL::Partition_is_valid_traits_2<Traits, PolygonIsValid>`

\cgalHeading{Implementation}

This test requires \f$ O(n)\f$ time for a polygon with \f$ n\f$ vertices.

*/
template< typename Traits >
class Is_y_monotone_2 {
public:

/// \name Creation
/// @{

/*!
`Traits` is a model of
the concept `IsYMonotoneTraits_2`
*/
Is_y_monotone_2(const Traits& t);

/// @}

/// \name Operations
/// @{

/*!

returns `true` iff the points of type `Traits::Point_2` in the
range [`first`,`beyond`) define a \f$ y\f$-monotone polygon.

*/
template<class InputIterator>
bool operator()(InputIterator first, InputIterator beyond);

/// @}

}; /* end Is_y_monotone_2 */
} /* end namespace CGAL */
