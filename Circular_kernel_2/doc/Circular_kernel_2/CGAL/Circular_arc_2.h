
namespace CGAL {

/*!
\ingroup PkgCircularKernel2GeometricClasses

\cgalModels{CircularKernel::CircularArc_2}

\sa `CGAL::Circular_arc_point_2<CircularKernel>`
\sa `CGAL::Line_arc_2<CircularKernel>`

*/
template< typename CircularKernel >
class Circular_arc_2 {
public:

/// \name Creation
/// @{

/*!
Constructs an arc from a full circle.
*/
Circular_arc_2(const Circle_2<CircularKernel> &c);

/*!
Constructs the circular arc supported by `c`, whose source is
`p` and whose target is `q` when traversing the circle in
counterclockwise direction.
\pre `p` and `q` lie on `c`.
*/
Circular_arc_2(const Circle_2<CircularKernel> &c,
const Circular_arc_point_2<CircularKernel> &p,
const Circular_arc_point_2<CircularKernel> &q);

/*!
Constructs an arc that is supported by the circle of type
`Circle_2<CircularKernel>` passing through the points `p`,
`q` and `r`. The source and target are respectively `p`
and `r`, when traversing the supporting circle in the
counterclockwise direction.
Note that, depending on the orientation of the point triple
`(p,q,r)`, `q` may not lie on the arc.
\pre `p`, `q`, and `r` are not collinear.
*/
Circular_arc_2(const Point_2<CircularKernel> &p,
const Point_2<CircularKernel> &q,
const Point_2<CircularKernel> &r);

/// @}

/// \name Access Functions
/// @{

/*!

*/
Circle_2<CircularKernel> supporting_circle();

/*!

returns the center of the supporting circle.
*/
Point_2<CircularKernel> const& center( ) const;

/*!

returns the squared radius of the supporting circle.
*/
CircularKernel::FT const& squared_radius( ) const;

/// @}

/// \name
/// A circular arc is not oriented. Still, its source and target endpoints can be defined, supposing that its supporting circle in traversed the counterclockwise direction from `source` to `target`.
/// @{

/*!

*/
Circular_arc_point_2<CircularKernel> source();

/*!

*/
Circular_arc_point_2<CircularKernel> target();

/// @}

/// \name
/// When the methods `source` and `target` return the same point, then
/// the arc is in fact a full circle.
///
/// When an arc is x-monotone, its left and right points can be
/// accessed directly:

/// @{

/*!
\pre `ca`.`is_x_monotone()`.
*/
Circular_arc_point_2<CircularKernel> left();

/*!
\pre `ca`.`is_x_monotone()`.
*/
Circular_arc_point_2<CircularKernel> right();

/*!
Returns a bounding box containing the arc.
*/
Bbox_2 bbox() const;

/// @}

/// \name Query Functions
/// @{

/*!
Tests whether the arc is x-monotone.
*/
bool is_x_monotone();

/*!
Tests whether the arc is y-monotone.
*/
bool is_y_monotone();
/// @}

}; /* end Circular_arc_2 */


/*!
Test for equality. Two arcs are equal, iff their non-oriented
supporting circles are equal (i.e.\ they have same center and same
squared radius) and their endpoints are equal.
\relates Circular_arc_2
*/
bool operator==(const Circular_arc_2<CircularKernel> &ca1, const Circular_arc_2<CircularKernel> &ca2);

/*!
Test for non-equality
\relates Circular_arc_2
*/
bool operator!=(const Circular_arc_2<CircularKernel> &ca1, const Circular_arc_2<CircularKernel> &ca2);

/*!

\relates Circular_arc_2
*/
istream& operator>> (std::istream& is, Circular_arc_2 & ca);


/*!

\relates Circular_arc_2
*/
ostream& operator<< (std::ostream& os, const Circular_arc_2 & ca);
} /* end namespace CGAL */
