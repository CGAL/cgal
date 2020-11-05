namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Classify a circle according to `sphere`, as defined in Section \ref sectionSKobjects.
\pre `c` lies on `sphere`.

\sa `CGAL::Circle_type`
*/

template <class SphericalKernel>
CGAL::Circle_type
classify(const CGAL::Circle_3<SphericalKernel> & c, const CGAL::Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Compares the \f$ \theta\f$-coordinates of `p` and `q` relatively to `sphere`.
\pre `p` and `q` lie on `sphere`, but do not coincide with the poles of `sphere`.

\sa \link compare_x_grp `CGAL::compare_x()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_x_at_y_grp `CGAL::compare_x_at_y()` \endlink
\sa \link compare_y_grp `CGAL::compare_y()` \endlink
\sa \link compare_yx_grp `CGAL::compare_yx()` \endlink
\sa \link compare_y_at_x_grp `CGAL::compare_y_at_x()` \endlink
\sa \link compare_z_grp `CGAL::compare_z()` \endlink
\sa `CGAL::compare_theta_z()`
*/
template <class SphericalKernel>
Comparison_result
compare_theta(const CGAL::Circular_arc_point_3<SphericalKernel> & p,
const CGAL::Circular_arc_point_3<SphericalKernel> & q,const CGAL::Sphere_3<SphericalKernel>& sphere);

/*!
\ingroup PkgCircularKernel3GeometricFunctions


Compares the \f$ \theta\f$-coordinates of `p` and of the meridian defined by `m` (see Section \ref sectionSKobjects)
in the cylindrical coordinate system relative to `sphere`.
\pre `p` lies on `sphere`, but does not coincide with its poles. `m` \f$ \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of `m` is \f$ 0\f$.

\sa \link compare_x_grp `CGAL::compare_x()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_x_at_y_grp `CGAL::compare_x_at_y()` \endlink
\sa \link compare_y_grp `CGAL::compare_y()` \endlink
\sa \link compare_yx_grp `CGAL::compare_yx()` \endlink
\sa \link compare_y_at_x_grp `CGAL::compare_y_at_x()` \endlink
\sa \link compare_z_grp `CGAL::compare_z()` \endlink
\sa `CGAL::compare_theta_z()`
*/
template <class SphericalKernel>
Comparison_result
compare_theta(const CGAL::Circular_arc_point_3<SphericalKernel> &p, const CGAL::Vector_3<SphericalKernel> &m, const CGAL::Sphere_3<SphericalKernel>& sphere );

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Compares the \f$ \theta\f$-coordinates of the meridian defined by `m` and  of `p` (see Section \ref sectionSKobjects)
in the cylindrical coordinate system relative to `sphere`.
\pre `p` lies on `sphere`, but does not coincide with its poles. `m` \f$ \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of `m` is \f$ 0\f$.

\sa \link compare_x_grp `CGAL::compare_x()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_x_at_y_grp `CGAL::compare_x_at_y()` \endlink
\sa \link compare_y_grp `CGAL::compare_y()` \endlink
\sa \link compare_yx_grp `CGAL::compare_yx()` \endlink
\sa \link compare_y_at_x_grp `CGAL::compare_y_at_x()` \endlink
\sa \link compare_z_grp `CGAL::compare_z()` \endlink
\sa `CGAL::compare_theta_z()`
*/
template <class SphericalKernel>
Comparison_result
compare_theta(const CGAL::Vector_3<SphericalKernel> &m,const CGAL::Circular_arc_point_3<SphericalKernel> &p, const CGAL::Sphere_3<SphericalKernel>& sphere );

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Compares `p` and `q` according to the lexicographic ordering on \f$ \theta\f$ and \f$ z\f$-coordinates
in the cylindrical coordinate system relative to `sphere`.
\pre `p` and `q` lie on `sphere`, but do not coincide with the poles of `sphere`.

\sa \link compare_x_grp `CGAL::compare_x()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_xy_grp `CGAL::compare_xy()` \endlink
\sa \link compare_x_at_y_grp `CGAL::compare_x_at_y()` \endlink
\sa \link compare_y_grp `CGAL::compare_y()` \endlink
\sa \link compare_yx_grp `CGAL::compare_yx()` \endlink
\sa \link compare_y_at_x_grp `CGAL::compare_y_at_x()` \endlink
\sa \link compare_z_grp `CGAL::compare_z()` \endlink
\sa `CGAL::compare_theta()`
*/
template <class SphericalKernel>
bool
compare_theta_z(const CGAL::Circular_arc_point_3<SphericalKernel> & p,
const CGAL::Circular_arc_point_3<SphericalKernel> & q, const CGAL::Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Tests whether the arc `a` is \f$ \theta\f$-monotone, i.e.\ the intersection of
any meridian anchored at the poles `sphere` and the arc `a`
is reduced to at most one point in general, and two points if a pole of `sphere` is an endpoint of the arc.
Note that a bipolar circle has no such arcs.
\pre `a` lies on `sphere`, and the supporting circle of `a` is not bipolar.

*/
template <class SphericalKernel>
bool
is_theta_monotone(const CGAL::Circular_arc_3<SphericalKernel> & a,const CGAL::Sphere_3<SphericalKernel>& sphere);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Returns the point on the circle that is extremal in \f$ \theta\f$ using the cylindrical coordinate system
relative to `sphere`, and that has the smallest (resp.\ largest)
\f$ \theta\f$-coordinate of the two points if `b` is `true` (resp.\ `false`).
See Section \ref sectionSKobjects for definitions.
\pre `c` lies on `sphere` and is a normal circle.

*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
theta_extremal_point(const CGAL::Circle_3<SphericalKernel> & c, const CGAL::Sphere_3<SphericalKernel> sphere, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ \theta\f$-extremal points of the
circle relatively to `sphere`. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted in the cylindrical coordinate system relative to `sphere`.
See Section \ref sectionSKobjects for definitions.
\pre `c` lies on `sphere` and is a normal circle.

*/
template < class SphericalKernel, class OutputIterator >
OutputIterator
theta_extremal_points(const CGAL::Circle_3<SphericalKernel> & c, const CGAL::Sphere_3<SphericalKernel>& sphere,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Returns the point on the sphere that is extremal in the
\f$ x\f$-direction, and that is the smallest (resp.\ largest) of the two
\f$ x\f$-extremal points for the lexicographic order if `b` is `true`
(resp.\ `false`).
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
x_extremal_point(const CGAL::Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ x\f$-axis.
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
x_extremal_point(const CGAL::Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ x\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template < class SphericalKernel, class OutputIterator >
OutputIterator
x_extremal_points(const CGAL::Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ x\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ x\f$-axis.
*/
template < class SphericalKernel, class OutputIterator >
OutputIterator
x_extremal_points(const CGAL::Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Returns the point on the sphere that is extremal in the
\f$ y\f$-direction, and that is the smallest (resp.\ largest) of the two
\f$ y\f$-extremal points for the lexicographic order if `b` is `true`
(resp.\ `false`).
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
y_extremal_point(const CGAL::Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ y\f$-axis.
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
y_extremal_point(const CGAL::Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ y\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template <class SphericalKernel class OutputIterator >
OutputIterator
y_extremal_points(const CGAL::Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ y\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ y\f$-axis.
*/
template < class SphericalKernel,class OutputIterator >
OutputIterator
y_extremal_points(const CGAL::Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Returns the point on the sphere that is extremal in the
\f$ z\f$-direction, and that is the smallest (resp.\ largest) of the two
\f$ z\f$-extremal points for the lexicographic order if `b` is `true`
(resp.\ `false`).
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
z_extremal_point(const CGAL::Sphere_3<SphericalKernel> & c, bool b);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Same for a circle.
\pre The circle is not contained in a plane orthogonal to the \f$ z\f$-axis.
*/
template <class SphericalKernel>
CGAL::Circular_arc_point_3<SphericalKernel>
z_extremal_point(const CGAL::Circle_3<SphericalKernel> & c, bool b);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ z\f$-extremal points of the
sphere. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically sorted.
*/
template < class SphericalKernel, class OutputIterator >
OutputIterator
z_extremal_points(const CGAL::Sphere_3<SphericalKernel> & c,
OutputIterator res);

/*!
\ingroup PkgCircularKernel3GeometricFunctions

Copies in the output iterator the \f$ z\f$-extremal points of the
circle. `res` iterates on elements of type
`Circular_arc_point_3<SphericalKernel>`, lexicographically
sorted.
\pre The circle is not contained in a plane orthogonal to the \f$ z\f$-axis.
*/
template < class SphericalKernel, class OutputIterator >
OutputIterator
z_extremal_points(const CGAL::Circle_3<SphericalKernel> & c,
OutputIterator res);

} /* namespace CGAL */

