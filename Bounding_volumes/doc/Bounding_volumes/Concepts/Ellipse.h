/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

An object `ellipse` of the class `Ellipse` is an ellipse in two-dimensional
Euclidean plane \f$ \E^2\f$. Its boundary splits the plane into a bounded
and an unbounded side. By definition, an empty `ellipse` has no
boundary and no bounded side, i.e.\ its unbounded side equals the
whole plane \f$ \E^2\f$.
*/
class Ellipse {
public:
/// \name Types
/// @{

/*!
Point type.
*/
typedef unspecified_type Point;

/// @}

/// \name Creation
/// @{

/*!
sets `ellipse` to the empty ellipse.
*/
void  set( );


/*!
sets `ellipse` to the ellipse containing exactly `p`.
*/
void  set( const Point& p);


/*!

sets `ellipse` to the ellipse containing exactly the segment
connecting `p` and `q`. The algorithm guarantees
that `set` is never called with two equal points.

*/
void  set( const Point& p,
                              const Point& q);


/*!

sets `ellipse` to the smallest ellipse through `p`,`q`,`r`.
The algorithm guarantees that `set` is never called with
three collinear points.
*/
void  set( const Point& p,
                              const Point& q,
                              const Point& r);


/*!

sets `ellipse` to the smallest ellipse through
`p`,`q`,`r`,`s`. The algorithm guarantees that
this ellipse exists.
*/
void  set( const Point& p,
                              const Point& q,
                              const Point& r,
                              const Point& s);


/*!

sets `ellipse` to the unique conic through
`p`,`q`,`r`,`s`,`t`. The algorithm
guarantees that this conic is an ellipse.
*/
void  set( const Point& p,
                              const Point& q,
                              const Point& r,
                              const Point& s,
                              const Point& t);

/// @}


/// \name Predicates
/// @{

/*!

returns `true`, iff `p` lies properly outside of `ellipse`.
*/
bool  has_on_unbounded_side( const Point& p) const;



/*!

returns `CGAL::ON_BOUNDED_SIDE`,
`CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies properly
inside, on the boundary, or properly outside of `ellipse`, resp.

\note Only needed, if the corresponding predicate of `CGAL::Min_ellipse_2` is used.
*/
CGAL::Bounded_side
                   bounded_side( const Point& p) const;


/*!

returns `true`, iff `p` lies properly inside `ellipse`.

\note Only needed, if the corresponding predicate of `CGAL::Min_ellipse_2` is used.
*/
bool  has_on_bounded_side( const Point& p) const;


/*!

returns `true`, iff `p` lies on the boundary
of `ellipse`.

\note Only needed, if the corresponding predicate of `CGAL::Min_ellipse_2` is used.
*/
bool  has_on_boundary( const Point& p) const;


/*!

returns `true`, iff `ellipse` is empty (this implies
degeneracy).

\note Only needed, if the corresponding predicate of `CGAL::Min_ellipse_2` is used.
*/
bool  is_empty( ) const;


/*!
returns `true`, iff `ellipse` is degenerate, i.e.\ if
`ellipse` is empty or equal to a single point.

\note Only needed, if the corresponding predicate of `CGAL::Min_ellipse_2` is used.
*/
bool  is_degenerate( ) const;

/// @}
/// \name I/O
/// The following I/O operator is only needed, if the corresponding I/O
/// operator of `CGAL::Min_ellipse_2` is used.
/// @{

/*!
writes `ellipse` to output stream `os`.

\note Only needed, if the corresponding I/O
operator of `CGAL::Min_ellipse_2` is used.

*/
ostream& operator<<(ostream& os, const Ellipse& ellipse);

/// @}

};
