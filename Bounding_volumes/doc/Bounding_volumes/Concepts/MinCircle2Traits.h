
/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

This concept defines the requirements for traits classes of
`CGAL::Min_circle_2<Traits>`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Min_circle_2_traits_2<K>}
\cgalHasModelsEnd

\sa `CGAL::Min_circle_2<Traits>`

*/

class MinCircle2Traits {
public:

/// \name Types
/// @{

/*!

The point type must provide default and copy constructor,
assignment and equality test.
*/
typedef unspecified_type Point;

/*!

The circle type must fulfill the requirements listed below
in the next section.
*/
typedef unspecified_type Circle;

/// @}

/// \name Variables
/// @{

/*!

The current circle. This variable is maintained by the algorithm,
the user should neither access nor modify it directly.
*/
Circle circle;

/// @}

/// \name Creation
/// Only default and copy constructor are required.
/// @{

/*!

*/
MinCircle2Traits( );

/*!

*/
MinCircle2Traits( const MinCircle2Traits&);

/// @}

/// \name Operations
/// The following predicate is only needed, if the member function
/// `is_valid` of `Min_circle_2` is used.
/// @{

/*!

returns constants `CGAL::LEFT_TURN`,
`CGAL::COLLINEAR`, or
`CGAL::RIGHT_TURN` iff `r` lies properly to the
left of, on, or properly to the right of the oriented line
through `p` and `q`, resp.
*/
CGAL::Orientation
orientation( const Point& p,
const Point& q,
const Point& r) const;

/// @}

}; /* end MinCircle2Traits */

