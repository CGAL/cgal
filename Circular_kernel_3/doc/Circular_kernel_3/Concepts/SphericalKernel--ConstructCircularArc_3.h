
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

The circular arc constructed from a circle, a source, and a target, is
defined as the set of points of the circle that lie between the source `p1`
and the target `p2`, when traversing the circle counterclockwise
seen from the side of the plane of the circle pointed by its <I>positive</I> normal
vectors.

In this definition, we say that a normal vector \f$ (a,b,c)\f$  is <I>positive</I> if
\f$ (a,b,c)>(0,0,0)\f$ (i.e.\ \f$ (a>0) || (a==0) \&\& (b>0) || (a==0)\&\&(b==0)\&\&(c>0)\f$).

*/

class SphericalKernel::ConstructCircularArc_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs an arc from a full circle.
*/
SphericalKernel::Circular_arc_3 operator()
(const SphericalKernel::Circle_3 &c);

/*!
Constructs the circular arc supported by `c`, whose source and target
are `p` and `q`, respectively.
\pre `p` and `q` lie on `c` and they are different.
*/
SphericalKernel::Circular_arc_3 operator()
(const SphericalKernel::Circle_3 &c,
const SphericalKernel::Circular_arc_point_3 &p,
const SphericalKernel::Circular_arc_point_3 &q);

/*!
Constructs an arc that is supported by the circle of type
`SphericalKernel::Circle_3` passing through the points `p`,
`q` and `r`. The source and target are respectively `p`
and `r`, when traversing the supporting circle in the
counterclockwise direction
seen from the side of the plane containing the circle pointed by its <I>positive</I>
normal vectors.
the circle.
Note that, depending on the orientation of the point triple
`(p,q,r)`, `q` may not lie on the arc.
\pre `p`, `q`, and `r` are not collinear.
*/
SphericalKernel::Circular_arc_3 operator()
(const SphericalKernel::Point_3 &p,
const SphericalKernel::Point_3 &q,
const SphericalKernel::Point_3 &r);

/// @}

}; /* end SphericalKernel::ConstructCircularArc_3 */

