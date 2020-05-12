namespace CGAL {

/*!
\ingroup PkgCircularKernel2GeometricFunctions


Checks whether the point lies in the vertical range defined by the
arc.
*/
bool
has_in_x_range(const Circular_arc_2<CircularKernel> &ca,
const Circular_arc_point_2<CircularKernel> &p);

/*!
\ingroup PkgCircularKernel2GeometricFunctions


Checks whether the point lies in the vertical range defined by the
line segment.
*/
bool
has_in_x_range(const Line_arc_2<CircularKernel> &ca,
const Circular_arc_point_2<CircularKernel> &p);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Checks whether the point lies on the circle.
*/
bool
has_on(const Circle_2<CircularKernel> &c,
const Circular_arc_point_2<CircularKernel> &p);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Copies in the output iterator the `x`-monotone sub-arcs of `ca`.
*/
template < class OutputIterator >
OutputIterator
make_x_monotone(const Circular_arc_2<CircularKernel> &ca,
OutputIterator res);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Copies in the output iterator the `xy`-monotone sub-arcs of `ca`.
*/
template < class OutputIterator >
OutputIterator
make_xy_monotone(const Circular_arc_2<CircularKernel> &ca,
OutputIterator res);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Returns the leftmost (resp.\ rightmost) point of the circle if `b` is
`true` (resp.\ `false`).
*/
Circular_arc_point_2<CircularKernel>
x_extremal_point(const Circle_2<CircularKernel> & c, bool b);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Copies in the output iterator the `x`-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_2<CircularKernel>`, sorted in `x`.
*/
template < class OutputIterator >
OutputIterator
x_extremal_points(const Circle_2<CircularKernel> & c,
OutputIterator res);


/*!
\ingroup PkgCircularKernel2GeometricFunctions


Returns the bottommost (resp.\ topmost) point of the circle if `b` is
`true` (resp.\ `false`).
*/
Circular_arc_point_2<CircularKernel>
y_extremal_point(const Circle_2<CircularKernel> & c, bool b);


/*!
\ingroup PkgCircularKernel2GeometricFunctions

Copies in the output iterator the `y`-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_2<CircularKernel>`, sorted in `y`.
*/
template < class OutputIterator >
OutputIterator
y_extremal_points(const Circle_2<CircularKernel> & c,
OutputIterator res);

/*!
\ingroup PkgCircularKernel2GeometricFunctions


Compares vertically the two arcs, to the right of the point `p`,
\pre `p` is an intersection point of the arcs, and the arcs are defined to the right of `p`.
*/
CGAL::Comparison_result
compare_y_to_right(const Circular_arc_2<CircularKernel> &ca1,
const Circular_arc_2<CircularKernel> &ca2,
Circular_arc_point_2<CircularKernel> &p);


} /* namespace CGAL */

