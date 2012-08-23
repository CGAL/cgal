namespace CGAL {

/*!
\ingroup PkgKernel23

An object of type `Circle_3` is a circle in the 
three-dimensional Euclidean space \f$ \E^3\f$. Note that the 
circle can be degenerate, i.e. the squared radius may be zero. 

\sa `Kernel::Circle_3` 

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
\pre `center` lies in `plane` and 		`sq_r` \f$ \geq\f$ 0. 
*/ 
Circle_3(Point_3<Kernel> const& center, 
Kernel::FT const& sq_r, 
Plane_3<Kernel> const& plane); 

/*! 
introduces a variable `c` of type `Circle_3`. 
It is initialized to the circle of center `center` and 
squared radius `sq_r` in a plane normal to 
the vector `n`. 
\pre `sq_r` \f$ \geq\f$ 0. 
*/ 
Circle_3(Point_3<Kernel> const& center, 
Kernel::FT const& sq_r, 
Vector_3<Kernel> const& n); 

/*! 
introduces a variable `c` of type `Circle_3`. 
It is initialized to the circle passing through the three points. 
\pre The three points are not collinear. 
*/ 
Circle_3(Point_3<Kernel> const& p, 
Point_3<Kernel> const& q, Point_3<Kernel> const& r); 

/*! 
introduces a variable `c` of type `Circle_3`. 
It is initialized to the circle along which the two spheres intersect. 
\pre The two spheres intersect along a circle. 
*/ 
Circle_3(Sphere_3<Kernel> const& sphere1, 
Sphere_3<Kernel> const& sphere2); 

/*! 
introduces a variable `c` of type `Circle_3`. 
It is initialized to the circle along which the sphere and the 
plane intersect. 
\pre The sphere and the plane intersect along a circle. 
*/ 
Circle_3(Sphere_3<Kernel> const& sphere, 
Plane_3<Kernel> const& plane); 

/*! 
introduces a variable `c` of type `Circle_3`. 
It is initialized to the circle along which the sphere and the 
plane intersect. 
\pre The sphere and the plane intersect along a circle. 
*/ 
Circle_3(Plane_3<Kernel> const& plane, 
Sphere_3<Kernel> const& sphere); 

/// @} 

/// \name Access Functions 
/// @{

/*! 

returns the center of `c`. 
*/ 
Point_3<Kernel> const& center( ) const; 

/*! 

returns the squared radius of `c`. 
*/ 
Kernel::FT const& squared_radius( ) const; 

/*! 

returns the supporting plane of `c`. 
*/ 
Plane_3<Kernel> const& supporting_plane( ) const; 

/*! 

returns the diametral sphere of `c`. 
*/ 
Sphere_3<Kernel> const& diametral_sphere( ) const; 

/*! 

returns the area of `c`, divided by \f$ \pi\f$. 
*/ 
Kernel::FT const& area_divided_by_pi( ) const; 

/*! 

returns an approximation of the area of `c`. 
*/ 
double const& approximate_area( ) const; 

/*! 

returns the squared length of `c`, divided by \f$ \pi^2\f$. 
*/ 
Kernel::FT const& squared_length_divided_by_pi_square( ) const; 

/*! 

returns an approximation of the squared length (i.e. perimeter) of `c`. 
*/ 
double const& approximate_squared_length( ) const; 

/// @} 

/// \name Predicates 
/// @{

/*! 

*/ 
bool has_on(Point_3<Kernel> const& p) const; 

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
i.e. if they have the same center, the same squared radius 
and the same supporting plane. 
\relates Circle_3 
*/ 
bool operator == (Circle_3<Kernel> const& c1, 
Circle_3<Kernel> const& c2); 

/*! 

\relates Circle_3 
*/ 
bool operator != (Circle_3<Kernel> const& c1, 
Circle_3<Kernel> const& c2); 

} /* end namespace CGAL */
