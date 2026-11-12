namespace CGAL {

/*!
\ingroup kernel_classes3

An object of type `Ball_3` is a ball in the
three-dimensional Euclidean space \f$ \E^3\f$. The ball is oriented, i.e.\ its
boundary has clockwise or counterclockwise orientation. The
boundary splits \f$ \E^3\f$ into a positive and a negative side, where the
positive side is to the left of the boundary. The boundary also
splits \f$ \E^3\f$ into a bounded and an unbounded side. Note that the
ball can be degenerated, i.e.\ the squared radius may be zero.

\cgalModels{Kernel::Ball_3,Hashable if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`}

*/
template< typename Kernel >
class Ball_3 {
public:

/// \name Creation
/// @{

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the ball with center `center`,
squared radius `squared_radius` and orientation
`orientation`.
\pre `orientation != COPLANAR` and `squared_radius >= 0`.

\cgalEpicExact
*/
Ball_3( const Point_3<Kernel> & center,
          const Kernel::FT & squared_radius,
          const Orientation & orientation = COUNTERCLOCKWISE);

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the unique ball which passes through
the points `p`, `q`, `r` and `s`. The orientation of
the ball is the orientation of the point quadruple `p`,
`q`, `r`, `s`.
\pre `p`, `q`, `r`, and `s` are not coplanar.
*/
Ball_3( const Point_3<Kernel> & p,
const Point_3<Kernel> & q,
const Point_3<Kernel> & r,
const Point_3<Kernel> & s);

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the smallest ball which passes through
the points `p`, `q`, and `r`. The orientation of
the ball is `o`. \pre `o != COPLANAR`.
*/
Ball_3( const Point_3<Kernel> & p,
const Point_3<Kernel> & q,
const Point_3<Kernel> & r,
const Orientation& o = COUNTERCLOCKWISE);

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the smallest ball which passes through
the points `p` and `q`. The orientation of
the ball is `o`. \pre `o != COPLANAR`.
*/
Ball_3( const Point_3<Kernel> & p,
const Point_3<Kernel> & q,
const Orientation& o = COUNTERCLOCKWISE);

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the ball with center `center`, squared
radius zero and orientation `orientation`.
\pre `orientation != COPLANAR`.
\post `c.is_degenerate()` = `true`.

\cgalEpicExact
*/
Ball_3( const Point_3<Kernel> & center,
          const Orientation& orientation = COUNTERCLOCKWISE);

/*!

introduces a variable `c` of type `Ball_3`.
It is initialized to the diametral ball of the circle.
\cgalEpicExact
*/
Ball_3( const Circle_3<Kernel> & c );

/// @}

/// \name Access Functions
/// @{

/*!

returns the center of `c`.
\cgalEpicExact
*/
const Point_3<Kernel> & center( ) const;

/*!

returns the squared radius of `c`.
\cgalEpicExact
*/
Kernel::FT const& squared_radius( ) const;

/*!

returns the orientation of `c`.
\cgalEpicExact
*/
Orientation const& orientation( ) const;

/*!

returns `true`, iff `c` and `ball2` are equal,
i.e.\ if they have the same center, same squared radius and
same orientation.
*/
bool operator == ( const Ball_3<Kernel> & ball2) const;

/*!

returns `true`, iff `c` and `ball2` are not equal.
*/
bool operator != (  const Ball_3<Kernel> & ball2) const;

/// @}

/// \name Predicates
/// @{

/*!

returns `true`, iff `c` is degenerate, i.e.\ if `c` has squared radius zero.
*/
bool is_degenerate( ) const;

/*!

returns either the constant \ref ON_ORIENTED_BOUNDARY,
\ref ON_POSITIVE_SIDE, or \ref ON_NEGATIVE_SIDE,
iff `p` lies on the boundary, properly on the
positive side, or properly on the negative side
of `c`, resp.
*/
Oriented_side
oriented_side( const Point_3<Kernel> & p) const;

/*!

returns \ref ON_BOUNDED_SIDE,
\ref ON_BOUNDARY, or \ref ON_UNBOUNDED_SIDE
iff `p` lies properly inside, on the boundary, or properly
outside of `c`, resp.
*/
Bounded_side
bounded_side( const Point_3<Kernel> & p) const;

/*!

*/
bool has_on_positive_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_negative_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_boundary(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_bounded_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_unbounded_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_unbounded_side(const Iso_cuboid_3<Kernel>& b) const

/*!

*/
bool has_on(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on(const Circle_3<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!

returns the ball with the same center and squared radius as
`c` but with opposite orientation.
\cgalEpicExact
*/
Ball_3<Kernel> opposite() const;

/*!

returns the ball obtained by applying `at` on `c`.
\pre `at` is an orthogonal transformation.
*/
Ball_3<Kernel> orthogonal_transform(
Aff_transformation_3<Kernel> const& at) const;

/*!

returns a bounding box containing `c`.
*/
Bbox_3 bbox() const;

/// @}

}; /* end Ball_3 */
} /* end namespace CGAL */
