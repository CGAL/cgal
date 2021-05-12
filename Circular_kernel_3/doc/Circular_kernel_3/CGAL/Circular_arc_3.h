
namespace CGAL {

/*!
\ingroup PkgCircularKernel3GeometricClasses

\cgalModels `SphericalKernel::CircularArc_3`

\cgalHeading{I/O}

The input/output format of a circular arc consists of the supporting circle
represented as a `Circle_3` object followed by the source and target
points of the arc represented as two `Circular_arc_point_3` objects.
The defined arc is the unique arc constructed from such three objects.

\sa `CGAL::Circular_arc_point_3<SphericalKernel>`
\sa `CGAL::Line_arc_3<SphericalKernel>`

*/
template< typename SphericalKernel >
class Circular_arc_3 {
public:

/// \name Creation
/// @{
/*!
Constructs an arc from a full circle.
*/
Circular_arc_3(const Circle_3<SphericalKernel> &c);

/*!
Constructs an arc from a full circle, using pt as source and target.
*/
Circular_arc_3(const Circle_3<SphericalKernel> &c, const Circular_arc_point_3& pt);

/// @}

/// The circular arc constructed from a circle, a source, and a
/// target, is defined as the set of points of the circle that lie
/// between the source `p1` and the target `p2`, when traversing the
/// circle counterclockwise seen from the side of the plane of the
/// circle pointed by its <I>positive</I> normal vectors.
///
/// In this
/// definition, we say that a normal vector \f$ (a,b,c)\f$ is <I>positive</I>
/// if \f$ (a,b,c)>(0,0,0)\f$ (i.e.\ \f$ (a>0) || (a==0) \&\& (b>0) || (a==0)\&\&(b==0)\&\&(c>0)\f$).
/// @{

/*!
Constructs the circular arc supported by `c`, whose source and target
are `p` and `q`, respectively.
\pre `p` and `q` lie on `c` and are different.
*/
Circular_arc_3(const Circle_3<SphericalKernel> &c,
const Circular_arc_point_3<SphericalKernel> &p,
const Circular_arc_point_3<SphericalKernel> &q);

/*!
Constructs an arc that is supported by the circle of type
`Circle_3<SphericalKernel>` passing through the points `p`,
`q` and `r`. The source and target are respectively `p`
and `r`, when traversing the supporting circle in the
counterclockwise direction
seen from the side of the plane containing the circle pointed by its <I>positive</I>
normal vectors.
Note that, depending on the orientation of the point triple
`(p,q,r)`, `q` may not lie on the arc.
\pre `p`, `q`, and `r` are not collinear.
*/
Circular_arc_3(const Point_3<SphericalKernel> &p,
const Point_3<SphericalKernel> &q,
const Point_3<SphericalKernel> &r);

/// @}

/// \name Access Functions
/// @{

/*!

*/
Circle_3<SphericalKernel> supporting_circle();

/*!

returns the center of the supporting circle.
*/
Point_3<SphericalKernel> const& center( ) const;

/*!

returns the squared radius of the supporting circle.
*/
SphericalKernel::FT const& squared_radius( ) const;

/*!

*/
Plane_3<SphericalKernel> supporting_plane();

/*!

*/
Sphere_3<SphericalKernel> diametral_sphere();

/// @}

/// \name
/// When the methods `source` and `target` return the same point, then
/// the arc is in fact a full circle. When the arc was constructed
/// from its (full) underlying circle, then source and target both
/// return the smallest \f$ x\f$-extremal point of the circle if the
/// circle is not contained in a plane \f$ x=A\f$, and the smallest
/// \f$ y\f$-extremal point otherwise.
/// @{

/*!

*/
Circular_arc_point_3<SphericalKernel> source();

/*!

*/
Circular_arc_point_3<SphericalKernel> target();

/*!
Test for equality. Two arcs are equal, iff their non-oriented
supporting planes are equal, and the centers and squared
radii of their respective supporting circles are equal, and their
sources and targets are equal.
*/
bool operator==(const Circular_arc_3<SphericalKernel> &a1, const Circular_arc_3<SphericalKernel> &a2);

/*!
Test for nonequality.
*/
bool operator!=(const Circular_arc_3<SphericalKernel> &a1,
const Circular_arc_3<SphericalKernel> &a2);

/// @}

}; /* end Circular_arc_3 */

/*!

\relates Circular_arc_3
*/
istream& operator>> (std::istream& is, Circular_arc_3 & ca);

/*!

\relates Circular_arc_3
*/
ostream& operator<< (std::ostream& os, const Circular_arc_3 & ca);


} /* end namespace CGAL */
