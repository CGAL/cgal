namespace CGAL {

/*!
\ingroup PkgCircularKernel2


Checks whether the point lies in the vertical range defined by the
arc.
*/
bool 
has_in_x_range(const Circular_arc_2<CircularKernel> &ca, 
const Circular_arc_point_2<CircularKernel> &p);

/*!
\ingroup PkgCircularKernel2


Checks whether the point lies in the vertical range defined by the
line segment.
*/
bool 
has_in_x_range(const Line_arc_2<CircularKernel> &ca, 
const Circular_arc_point_2<CircularKernel> &p);


/*!
\ingroup PkgCircularKernel2


Checks whether the point lies on the circle.
*/
bool
has_on(const Circle_2<CircularKernel> &c, 
const Circular_arc_point_2<CircularKernel> &p);


/*!
\ingroup PkgCircularKernel2


Copies in the output iterator the \f$ x\f$-monotone sub-arcs of \f$ ca\f$.
*/
template < class OutputIterator >
OutputIterator
make_x_monotone(const Circular_arc_2<CircularKernel> &ca, 
OutputIterator res);


/*!
\ingroup PkgCircularKernel2


Copies in the output iterator the \f$ xy\f$-monotone sub-arcs of \f$ ca\f$.
*/
template < class OutputIterator >
OutputIterator
make_xy_monotone(const Circular_arc_2<CircularKernel> &ca, 
OutputIterator res);


/*!
\ingroup PkgCircularKernel2


Returns the leftmost (resp. rightmost) point of the circle if \f$ b\f$ is
`true` (resp. `false`).
*/
Circular_arc_point_2<CircularKernel>
x_extremal_point(const Circle_2<CircularKernel> & c, bool b);


/*!
\ingroup PkgCircularKernel2


Copies in the output iterator the \f$ x\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_2<CircularKernel>`, sorted in \f$ x\f$.
*/
template < class OutputIterator >
OutputIterator
x_extremal_points(const Circle_2<CircularKernel> & c,
OutputIterator res);


/*!
\ingroup PkgCircularKernel2


Returns the bottommost (ressp. topmost) point of the circle if \f$ b\f$ is
`true` (resp. `false`).
*/
Circular_arc_point_2<CircularKernel>
y_extremal_point(const Circle_2<CircularKernel> & c, bool b);


/*!
\ingroup PkgCircularKernel2

Copies in the output iterator the \f$ y\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_2<CircularKernel>`, sorted in \f$ y\f$.
*/
template < class OutputIterator >
OutputIterator
y_extremal_points(const Circle_2<CircularKernel> & c,
OutputIterator res);

/*!
\ingroup PkgCircularKernel2


Compares vertically the two arcs, to the right of the point \f$ p\f$,
\pre \f$ p\f$ is an intersection point of the arcs, and the arcs are defined to the right of \f$ p\f$.
*/
CGAL::Comparison_result 
compare_y_to_right(const Circular_arc_2<CircularKernel> &ca1,
const Circular_arc_2<CircularKernel> &ca2,
Circular_arc_point_2<CircularKernel> &p);


} /* namespace CGAL */

