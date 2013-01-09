namespace CGAL {

/*!
\ingroup kernel_enums

converts between the various enums provided by the \cgal kernel.
The conversion preserves the order of the values.

\sa `CGAL::Sign` 
\sa `CGAL::Comparison_result` 
\sa `CGAL::Orientation` 
\sa `CGAL::Oriented_side` 
\sa `CGAL::Bounded_side` 
\sa `CGAL::Angle` 
\sa `CGAL::Uncertain<T>` 
*/
template < typename T, typename U >
T enum_cast(const U&u);

/*!
\ingroup kernel_enums

returns the opposite side (for example `CGAL::ON_POSITIVE_SIDE` if
`o`==`CGAL::ON_NEGATIVE_SIDE`), or `CGAL::ON_ORIENTED_BOUNDARY` if
`o`==`CGAL::ON_ORIENTED_BOUNDARY`.
*/
Oriented_side opposite(const Oriented_side &o);

/*!
\ingroup kernel_enums

returns the opposite side (for example `CGAL::ON_BOUNDED_SIDE` if
`o`==`CGAL::ON_UNBOUNDED_SIDE`), or returns `CGAL::ON_BOUNDARY` if
`o`==`CGAL::ON_BOUNDARY`.
*/
Bounded_side opposite(const Bounded_side &o);

/*!
\ingroup kernel_enums
\sa `angle_grp`
 */
enum Angle {OBTUSE, RIGHT, ACUTE};

/*!
\ingroup kernel_enums
\sa `CGAL::opposite(const Bounded_side& o)`
*/
enum Bounded_side {ON_UNBOUNDED_SIDE, ON_BOUNDARY, ON_BOUNDED_SIDE};

/*!
\ingroup kernel_enums
*/
enum Comparison_result { SMALLER, EQUAL, LARGER };

/*!
\ingroup kernel_enums
\sa `CGAL::Orientation`
*/
enum  Sign { NEGATIVE, ZERO, POSITIVE };

/*!
\ingroup kernel_enums
\sa `CGAL::LEFT_TURN` 
\sa `CGAL::RIGHT_TURN` 
\sa `CGAL::COLLINEAR`
\sa `CGAL::CLOCKWISE` 
\sa `CGAL::COUNTERCLOCKWISE`
\sa `CGAL::COPLANAR` 
*/
typedef Sign Orientation;

/*!
\ingroup kernel_enums
*/
enum Oriented_side {ON_NEGATIVE_SIDE, ON_ORIENTED_BOUNDARY, ON_POSITIVE_SIDE };

/*!
\ingroup kernel_enums

\sa `CGAL::COUNTERCLOCKWISE`
*/
const CGAL::Orientation CLOCKWISE = NEGATIVE;

/*!
\ingroup kernel_enums
\sa `CGAL::CLOCKWISE`
*/
const CGAL::Orientation COUNTERCLOCKWISE = POSITIVE;

/*!
\ingroup kernel_enums
\sa `CGAL::LEFT_TURN` 
\sa `CGAL::RIGHT_TURN` 
*/
const CGAL::Orientation COLLINEAR = ZERO;

/*!
\ingroup kernel_enums

\sa `CGAL::COLLINEAR` 

\sa `CGAL::RIGHT_TURN` 
*/
const CGAL::Orientation LEFT_TURN = POSITIVE;

/*!
\ingroup kernel_enums

\sa `CGAL::COLLINEAR`
\sa `CGAL::LEFT_TURN`

*/
const CGAL::Orientation RIGHT_TURN = NEGATIVE;

/*!
\ingroup kernel_enums
*/
const CGAL::Orientation COPLANAR  = ZERO;

/*!
\ingroup kernel_enums
*/
const CGAL::Orientation DEGENERATE = ZERO;

/*!
\ingroup kernel_enums

A symbolic constant used to construct zero length vectors.

\sa `CGAL::Vector_2<Kernel>` 
\sa `CGAL::Vector_3<Kernel>` 

*/
const CGAL::Null_vector NULL_VECTOR;

/*!
\ingroup kernel_enums

A symbolic constant which denotes the point at the origin.
This constant is used in the conversion between points and vectors.

\cgalHeading{Example}

\code
  Point_2< Cartesian<Exact_NT> >  p(1.0, 1.0), q;
  Vector2< Cartesian<Exact_NT> >  v;
  v = p - ORIGIN;
  q = ORIGIN + v;  
  assert( p == q );
\endcode

\sa `CGAL::Point_2<Kernel>` 
\sa `CGAL::Point_3<Kernel>` 

*/
const CGAL::Origin ORIGIN;

} /* namespace CGAL */
