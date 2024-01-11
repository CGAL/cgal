namespace CGAL {

/*!
\ingroup kernel_classes3

An object `c` of type `Circle_3` is a circle in the
three-dimensional Euclidean space \f$ \E^3\f$. Note that the
circle can be degenerate, i.e.\ the squared radius may be zero.

\cgalModels{Kernel::Circle_3}

*/
template< typename Kernel >
class Circle_3 {
public:

/// \name Creation
/// @{

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle of center `center` and
squared radius `sq_r` in plane `plane`.
\pre `center` lies in `plane` and `sq_r >= 0`.

\cgalEpicExact
*/
Circle_3(const Point_3<Kernel> &center,
         const Kernel::FT &sq_r,
         const Plane_3<Kernel> & plane);

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle of center `center` and
squared radius `sq_r` in a plane normal to
the vector `n`.
\pre `sq_r >= 0`.
*/
Circle_3(const Point_3<Kernel> & center,
         const Kernel::FT & sq_r,
         const Vector_3<Kernel> & n);

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle passing through the three points.
\pre The three points are not collinear.
*/
Circle_3(const Point_3<Kernel> & p,
         const Point_3<Kernel> & q, const Point_3<Kernel> & r);

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle along which the two spheres intersect.
\pre The two spheres intersect along a circle.
*/
Circle_3(const Sphere_3<Kernel> & sphere1,
         const Sphere_3<Kernel> & sphere2);

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle along which the sphere and the
plane intersect.
\pre The sphere and the plane intersect along a circle.
*/
  Circle_3(const Sphere_3<Kernel> & sphere,
           const Plane_3<Kernel> & plane);

/*!
introduces a variable `c` of type `Circle_3`.
It is initialized to the circle along which the sphere and the
plane intersect.
\pre The sphere and the plane intersect along a circle.
*/
Circle_3(const Plane_3<Kernel> & plane,
         const Sphere_3<Kernel> & sphere);

/// @}

/// \name Access Functions
/// @{

/*!

returns the center of `c`.
\cgalEpicExact
*/
const Point_3<Kernel> & center( ) const;

/*!

returns the squared radius of `c`.
\cgalEpicExact
*/
const Kernel::FT & squared_radius( ) const;

/*!

returns the supporting plane of `c`.
\cgalEpicExact
*/
const Plane_3<Kernel> & supporting_plane( ) const;

/*!

returns the diametral sphere of `c`.
\cgalEpicExact
*/
const Sphere_3<Kernel> & diametral_sphere( ) const;

/*!

returns the area of `c`, divided by \f$ \pi\f$.
*/
Kernel::FT const& area_divided_by_pi( ) const;

/*!

returns an approximation of the area of `c`.
*/
double approximate_area( ) const;

/*!

returns the squared length of `c`, divided by \f$ \pi^2\f$.
*/
Kernel::FT squared_length_divided_by_pi_square( ) const;

/*!

returns an approximation of the squared length (i.e.\ perimeter) of `c`.
*/
double approximate_squared_length( ) const;

/// @}

/// \name Predicates
/// @{

/*!

*/
bool has_on(const Point_3<Kernel> & p) const;

/// @}

/// \name Operations
/// @{

/*!

returns a bounding box containing `c`.
*/
Bbox_3 bbox() const;

/// @}

}; /* end Circle_3 */

/*!
returns `true`, iff `c1` and `c2` are equal,
i.e.\ if they have the same center, the same squared radius
and the same supporting plane.
\relates Circle_3
*/
bool operator == (const Circle_3<Kernel>& c1,
Circle_3<Kernel> const& c2);

/*!

\relates Circle_3
*/
bool operator != (const Circle_3<Kernel> & c1,
Circle_3<Kernel> const& c2);

} /* end namespace CGAL */
