/*!
\ingroup PkgGeneratorsConcepts
\cgalConcept

The concept `RandomConvexHullTraits_2` describes the requirements for the traits
class used by the function `random_convex_hull_in_disc_2()`.

\cgalHasModel \cgal kernels.

\cgalHeading{Operations}

The following two member functions returning instances of the above predicate
object types are required.

*/

class RandomConvexHullTraits_2 {
public:

/// \name Types
/// @{

/*!
The coordinate type of the points of the polygon.
*/
typedef unspecified_type FT;

/*!
The point type of the polygon.
*/
typedef unspecified_type Point_2;

/*!
The segment type of the polygon.
*/
typedef unspecified_type Segment_2;

/*!
Predicate object type that determines the orientation of three points.
It must provide `Orientation operator()(Point_2 p, Point_2 q,
Point_2 r)` that
returns `LEFT_TURN`, if \f$ r\f$ lies to the left of the oriented
line \f$ l\f$ defined by \f$ p\f$ and \f$ q\f$, returns `RIGHT_TURN` if \f$ r\f$
lies to the right of \f$ l\f$, and returns `COLLINEAR` if \f$ r\f$ lies
on \f$ l\f$.
*/
typedef unspecified_type Orientation_2;
/*!
A function object to compare the x-coordinate of two points.

Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
`SMALLER, EQUAL` or `LARGER`
according to the
\f$ x\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
A function object to compare the y-coordinate of two points.
Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
(`SMALLER, EQUAL` or `LARGER`)
according to the
\f$ y\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;


/// @}

/// \name Operations
/// @{

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();
/*!

*/
Orientation_2 orientation_2_object();


/// @}

}; /* end RandomConvexHullTraits_2 */
