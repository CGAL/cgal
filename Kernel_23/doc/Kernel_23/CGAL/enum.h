namespace CGAL {

/*!
\ingroup PkgKernel23

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
\ingroup PkgKernel23

returns the opposite side (for example `ON_POSITIVE_SIDE` if
`o==ON_NEGATIVE_SIDE`), or `ON_ORIENTED_BOUNDARY` if
`o==ON_ORIENTED_BOUNDARY`.
*/
Oriented_side opposite(const Oriented_side &o);

/*!
\ingroup PkgKernel23

returns the opposite side (for example `BOUNDED_SIDE` if
`o==UNBOUNDED_SIDE`), or returns `ON_BOUNDARY` if
`o==ON_BOUNDARY`.
*/
Bounded_side opposite(const Bounded_side &o);

/*!
\ingroup PkgKernel23
\sa `CGAL::angle`
 */
enum Angle {OBTUSE, RIGHT, ACUTE};

/*!
\ingroup PkgKernel23
\sa `CGAL::opposite(const Bounded_side& o)`
*/
enum Bounded_side {ON_UNBOUNDED_SIDE, ON_BOUNDARY, ON_BOUNDED_SIDE};

/*!
\ingroup PkgKernel23
*/
enum Comparison_result { SMALLER, EQUAL, LARGER };

/*!
\ingroup PkgKernel23
\sa CGAL::Orientation
*/
enum  Sign { NEGATIVE, ZERO, POSITIVE };

/*!
\ingroup PkgKernel23
\sa `CGAL::LEFT_TURN` 
\sa `CGAL::RIGHT_TURN` 
\sa `CGAL::COLLINEAR`
\sa `CGAL::CLOCKWISE` 
\sa `CGAL::COUNTERCLOCKWISE`
\sa `CGAL::COPLANAR` 
*/
typedef Sign Orientation;

/*!
\ingroup PkgKernel23
*/
enum Oriented_side {ON_NEGATIVE_SIDE, ON_ORIENTED_BOUNDARY, ON_POSITIVE_SIDE };

/*!
\ingroup PkgKernel23

\sa `CGAL::COUNTERCLOCKWISE`
*/
const Orientation CLOCKWISE = NEGATIVE;

/*!
\ingroup PkgKernel23
\sa `CGAL::CLOCKWISE`
*/
const Orientation COUNTERCLOCKWISE = POSITIVE;

/*!
\ingroup PkgKernel23
\sa `CGAL::LEFT_TURN` 
\sa `CGAL::RIGHT_TURN` 
*/
const Orientation COLLINEAR = ZERO;

/*!
\ingroup PkgKernel23

\sa `CGAL::COLLINEAR` 

\sa `CGAL::RIGHT_TURN` 
*/
const Orientation LEFT_TURN = POSITIVE;

/*!
\ingroup PkgKernel23

\sa `CGAL::COLLINEAR`
\sa `CGAL::LEFT_TURN`

*/
const Orientation RIGHT_TURN = NEGATIVE;

/*!
\ingroup PkgKernel23
*/
const Orientation COPLANAR  = ZERO;

/*!
\ingroup PkgKernel23
*/
const Orientation DEGENERATE = ZERO;

/*!
\ingroup PkgKernel23

A symbolic constant used to construct zero length vectors.

\sa `CGAL::Vector_2<Kernel>` 
\sa `CGAL::Vector_3<Kernel>` 

*/
const Null_vector NULL_VECTOR;

/*!
\ingroup PkgKernel23

A symbolic constant which denotes the point at the origin.
This constant is used in the conversion between points and vectors.

### Example ###

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
const Origin ORIGIN;

} /* namespace CGAL */
