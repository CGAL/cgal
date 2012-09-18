namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Classify a circle according to `sphere`, as defined in section \ref sectionSKobjects.
\pre `c` lies on `sphere`.

\sa `CGAL::Circle_type`
*/

CGAL::Circle_type
classify(const Circle_3<SphericalKernel> & c, const Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Compares the \f$ \theta\f$-coordinates of \f$ p\f$ and \f$ q\f$ relatively to `sphere`.
\pre `p` and `q` lie on `sphere`, but do not coincide with the poles of `sphere`.

\sa `CGAL::compare_x`
\sa `CGAL::compare_xy`
\sa `CGAL::compare_xyz`
\sa `CGAL::compare_x_at_y`
\sa `CGAL::compare_y`
\sa `CGAL::compare_yx`
\sa `CGAL::compare_y_at_x`
\sa `CGAL::compare_z`
\sa `CGAL::compare_theta_z`
*/
Comparison_result
compare_theta(const Circular_arc_point_3<SphericalKernel> & p,
const Circular_arc_point_3<SphericalKernel> & q,const Sphere_3<SphericalKernel>& sphere);

/*!
\ingroup PkgSphericalKernel3


Compares the \f$ \theta\f$-coordinates of \f$ p\f$ and of the meridian defined by \f$ m\f$ (see section \ref sectionSKobjects)
in the cylindrical coordinate system relative to `sphere` .
\pre `p` lies on `sphere`, but does not coincide with its poles. \f$ m \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of \f$ m\f$ is \f$ 0\f$.

\sa `CGAL::compare_x`
\sa `CGAL::compare_xy`
\sa `CGAL::compare_xyz`
\sa `CGAL::compare_x_at_y`
\sa `CGAL::compare_y`
\sa `CGAL::compare_yx`
\sa `CGAL::compare_y_at_x`
\sa `CGAL::compare_z`
\sa `CGAL::compare_theta_z`
*/
Comparison_result
compare_theta(const SphericalKernel::Circular_arc_point_3 &p, const SphericalKernel::Vector_3 &m, const SphericalKernel::Sphere_3& sphere );

/*!
\ingroup PkgSphericalKernel3

Same as previous, with opposite result.

\sa `CGAL::compare_x`
\sa `CGAL::compare_xy`
\sa `CGAL::compare_xyz`
\sa `CGAL::compare_x_at_y`
\sa `CGAL::compare_y`
\sa `CGAL::compare_yx`
\sa `CGAL::compare_y_at_x`
\sa `CGAL::compare_z`
\sa `CGAL::compare_theta_z`
*/
Comparison_result
compare_theta(const SphericalKernel::Vector_3 &m,const SphericalKernel::Circular_arc_point_3 &p);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Compares \f$ p\f$ and \f$ q\f$ according to the lexicographic ordering on \f$ \theta\f$ and \f$ z\f$-coordinates
in the cylindrical coordinate system relative to `sphere`.
\pre `p` and `q` lie on `sphere`, but do not coincide with the poles of `sphere`.

\sa `CGAL::compare_x`
\sa `CGAL::compare_xy`
\sa `CGAL::compare_xyz`
\sa `CGAL::compare_x_at_y`
\sa `CGAL::compare_y`
\sa `CGAL::compare_yx`
\sa `CGAL::compare_y_at_x`
\sa `CGAL::compare_z`
\sa `CGAL::compare_theta`
*/
bool
compare_theta_z(const Circular_arc_point_3<SphericalKernel> & p,
const Circular_arc_point_3<SphericalKernel> & q,const Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Tests whether the arc \f$ a\f$ is \f$ \theta\f$-monotone, i.e. the intersection of
any meridian anchored at the poles `sphere` and the arc \f$ a\f$
is reduced to at most one point in general, and two points if a pole of `sphere` is an endpoint of the arc.
Note that a bipolar circle has no such arcs.
\pre `a` lies on `sphere`, and the supporting circle of `a` is not bipolar.

*/
bool
is_theta_monotone(const Circular_arc_3<SphericalKernel> & a,const Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Returns the point on the circle that is extremal in \f$ \theta\f$ using the cylindrical coordinate system
relative to `sphere`, and that has the smallest (resp. largest)
\f$ \theta\f$-coordinate of the two points if \f$ b\f$ is `true` (resp. `false`).
See section \ref sectionSKobjects for definitions.
\pre `c` lies on `sphere` and is a normal circle.

*/
Circular_arc_point_3<SphericalKernel>
theta_extremal_point(const Circle_3<SphericalKernel> & c, const Sphere_3<SphericalKernel> sphere, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ \theta\f$-extremal points of the
circle relatively to `sphere`. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted in the cylindrical coordinate system relative to `sphere`.
See section \ref sectionSKobjects for definitions.
\pre `c` lies on `sphere` and is a normal circle.

*/
template < class OutputIterator >
OutputIterator
theta_extremal_points(const Circle_3<SphericalKernel> & c, const Sphere_3<SphericalKernel>& sphere,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Returns the point on the sphere that is extremal in the
\f$ x\f$-direction, and that is the smallest (resp. largest) of the two
\f$ x\f$-extremal points for the lexicographic order if \f$ b\f$ is `true` 
(resp. `false`).
*/
Circular_arc_point_3<SphericalKernel>
x_extremal_point(const Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgSphericalKernel3

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ x\f$-axis.
*/
Circular_arc_point_3<SphericalKernel>
x_extremal_point(const Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ x\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template < class OutputIterator >
OutputIterator
x_extremal_points(const Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ x\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ x\f$-axis.
*/
template < class OutputIterator >
OutputIterator
x_extremal_points(const Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Returns the point on the sphere that is extremal in the
\f$ y\f$-direction, and that is the smallest (resp. largest) of the two
\f$ y\f$-extremal points for the lexicographic order if \f$ b\f$ is `true` 
(resp. `false`).
*/
Circular_arc_point_3<SphericalKernel>
y_extremal_point(const Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgSphericalKernel3

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ y\f$-axis.
*/
Circular_arc_point_3<SphericalKernel>
y_extremal_point(const Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ y\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template < class OutputIterator >
OutputIterator
y_extremal_points(const Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ y\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ y\f$-axis.
*/
template < class OutputIterator >
OutputIterator
y_extremal_points(const Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Returns the point on the sphere that is extremal in the
\f$ z\f$-direction, and that is the smallest (resp. largest) of the two
\f$ z\f$-extremal points for the lexicographic order if \f$ b\f$ is `true` 
(resp. `false`).
*/
Circular_arc_point_3<SphericalKernel>
z_extremal_point(const Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgSphericalKernel3

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ z\f$-axis.
*/
Circular_arc_point_3<SphericalKernel>
z_extremal_point(const Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ z\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template < class OutputIterator >
OutputIterator
z_extremal_points(const Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgSphericalKernel3

Copies in the output iterator the \f$ z\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ z\f$-axis.
*/
template < class OutputIterator >
OutputIterator
z_extremal_points(const Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

