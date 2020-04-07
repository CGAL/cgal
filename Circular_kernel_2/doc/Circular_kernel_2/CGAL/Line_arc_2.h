
namespace CGAL {

/*!
\ingroup PkgCircularKernel2GeometricClasses

\cgalModels `CircularKernel::LineArc_2`

\cgalHeading{I/O}

The format for input/output is, for each line arc: a `Line_2`
(the supporting line) and two `Circular_arc_point_2` (the two endpoints),
under the condition that the endpoints are actually lying on the line.

\sa `CGAL::Circular_arc_point_2<CircularKernel>`
\sa `CGAL::Circular_arc_2<CircularKernel>`

*/
template< typename CircularKernel >
class Line_arc_2 {
public:

/// \name Creation
/// @{

/*!
Construct the line segment supported by `l`, whose source
is `p1` and whose target is `p2`.
\pre `p1` and `p2` lie on `l`.
*/
Line_arc_2(const Line_2<CircularKernel> &l,
const Circular_arc_point_2<CircularKernel> &p1,
const Circular_arc_point_2<CircularKernel> &p2);

/*!
Same.
*/
Line_arc_2(const Line_2<CircularKernel> &l,
const Point_2<CircularKernel> &p1,
const Point_2<CircularKernel> &p2);

/*!

*/
Line_arc_2(const Segment_2<CircularKernel> &s);

/// @}

/// \name Access Functions
/// @{

/*!

*/
Line_2<CircularKernel> supporting_line();

/*!

*/
Circular_arc_point_2<CircularKernel> source();

/*!

*/
Circular_arc_point_2<CircularKernel> target();

/*!

*/
Circular_arc_point_2<CircularKernel> left();

/*!

*/
Circular_arc_point_2<CircularKernel> right();

/*!
Returns a bounding box containing the line segment.
*/
Bbox_2 bbox() const;

/// @}

/// \name Query Functions
/// @{

/*!

*/
bool is_vertical();

/// @}

}; /* end Line_arc_2 */


/*!
Test for equality. Two arcs are equal, iff their non-oriented
supporting lines are equal (i.e.\ they contain the same set of
points) and their endpoints are equal.
\relates Line_arc_2
*/
bool operator==(const Line_arc_2<CircularKernel> &la1, const Line_arc_2<CircularKernel> &la2);

/*!
Test for non-equality.
\relates Line_arc_2
*/
bool operator!=(const Line_arc_2<CircularKernel> &la1, const Line_arc_2<CircularKernel> &la2);

/*!

\relates Line_arc_2
*/
istream& operator>> (std::istream& is, Line_arc_2 & ca);

/*!

\relates Line_arc_2
*/
ostream& operator<< (std::ostream& os, const Line_arc_2 & ca);

} /* end namespace CGAL */
