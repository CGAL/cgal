/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

An object of the class `Circle` is a circle in two-dimensional
Euclidean plane \f$ \E^2\f$. Its boundary splits the plane into a bounded
and an unbounded side. By definition, an empty `Circle` has no
boundary and no bounded side, i.e.\ its unbounded side equals the
whole plane \f$ \E^2\f$. A `Circle` containing exactly one point \f$ p\f$
has no bounded side, its boundary is \f$ \{p\}\f$, and its unbounded side
equals \f$ \E^2 \setminus \{p\}\f$.

*/
class Circle
{
public:
/// \name Types
/// @{
/*!
Point type.
*/
typedef unspecified_type Point;

/*!
Distance type. The function `squared_radius()` (see below)
returns an object of this type.

\note Only needed, if the member function `is_valid()`
of `Min_circle_2` is used.

*/
typedef unspecified_type Distance;

/// @}

/// \name Creation
/// @{

/**
 * sets `circle` to the empty circle.
 */
void set();

/**
 * sets `circle` to the circle containing exactly \f$ \{\mbox{`p` }\}\f$.
 */
void  set( const Point& p);


/*!
sets `circle` to the circle with diameter equal to the segment
connecting `p` and `q`. The algorithm guarantees that `set` is never
called with two equal points.
*/
void  set( const Point& p, const Point& q);

/*!
sets `circle` to the circle through `p`,`q`,`r`.  The algorithm
guarantees that `set` is never called with three collinear points.
*/
void  set( const Point& p, const Point& q, const Point& r);

/// @}

/// \name Predicates
/// @{

/*!
returns `true`, iff `p` lies properly outside of `circle`.
*/
bool  has_on_unbounded_side( const Point& p) const;

/*!

returns `CGAL::ON_BOUNDED_SIDE`, `CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies properly inside, on the
boundary, or properly outside of `circle`, resp.

\note Only needed, if the corresponding predicate of `Min_circle_2` is used.
*/
CGAL::Bounded_side
                   bounded_side( const Point& p) const;


/*!
returns `true`, iff `p` lies properly inside `circle`.

\note Only needed, if the corresponding predicate of `Min_circle_2` is used.
*/
bool  has_on_bounded_side( const Point& p) const;


/*!

returns `true`, iff `p` lies on the boundary of `circle`.

\note Only needed, if the corresponding predicate of `Min_circle_2` is used.
*/
bool  has_on_boundary( const Point& p) const;


/*!

returns `true`, iff `circle` is empty (this implies degeneracy).

\note Only needed, if the corresponding predicate of `Min_circle_2` is used.
*/
bool  is_empty( ) const;


/*!

returns `true`, iff `circle` is degenerate, i.e.\ if `circle` is empty
or equal to a single point.

\note Only needed, if the corresponding predicate of `Min_circle_2` is used.
*/
bool  is_degenerate( ) const;

/// @}

/// \name Additional Operations for Checking
/// @{



/*!

returns `true`, iff `circle` and `circle2` are equal.

\note Only needed, if the member function `is_valid()` of `Min_circle_2` is used.
*/
bool  operator == ( const Circle& circle2) const;


/*!
returns the center of `circle`.

\note Only needed, if the member function `is_valid()` of `Min_circle_2` is used.
*/
Point  center( ) const;


/*!

returns the squared radius of `circle`.

\note Only needed, if the member function `is_valid()` of `Min_circle_2` is used.
*/
Distance  squared_radius( ) const;

/// @}

/// \name I/O
/// @{
/*!
writes `circle` to output stream `os`.
*/
std::ostream& operator<<( std::ostream& os, const Circle& circle);

/// @}

};
