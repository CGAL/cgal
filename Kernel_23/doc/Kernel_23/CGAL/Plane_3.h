namespace CGAL {
/*!
\ingroup kernel_classes3

An object `h` of the data type `Plane_3` is an oriented 
plane in the three-dimensional Euclidean space \f$ \E^3\f$. It is defined 
by the set of points with %Cartesian coordinates \f$ (x,y,z)\f$ that satisfy 
the plane equation 

\f[h :\;  a\, x +b\, y +c\, z + d = 0.\f]

The plane splits \f$ \E^3\f$ in a <I>positive</I> and a <I>negative side</I>. 
A point `p` with %Cartesian coordinates \f$ (px, py, pz)\f$ is on the 
positive side of `h`, iff \f$ a\, px +b\, py +c\, pz + d > 0\f$. 
It is on the negative side, iff \f$ a\, px +b\, py\, +c\, pz + d < 0\f$. 

\sa `Kernel::Plane_3` 

*/
template< typename Kernel >
class Plane_3 {
public:

/// \name Creation 
/// @{

/*!
creates a plane `h` defined by the equation 
\f$ a\, px +b\, py +c\, pz + d = 0\f$. 
Notice that `h` is degenerate if 
\f$ a = b = c = 0\f$. 
*/ 
Plane_3(const Kernel::RT &a, 
const Kernel::RT &b, 
const Kernel::RT &c, 
const Kernel::RT &d); 

/*!
creates a plane `h` passing through the points `p`, 
`q` and `r`. The plane is oriented such that `p`, 
`q` and `r` are oriented in a positive sense 
(that is counterclockwise) when seen from the positive side of `h`. 
Notice that `h` is degenerate if the points are collinear. 
*/ 
Plane_3(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q, 
const Point_3<Kernel> &r); 

/*!
introduces a plane `h` that passes through point `p` and 
that is orthogonal to `v`. 
*/ 
Plane_3(const Point_3<Kernel> &p, 
const Vector_3<Kernel> &v); 

/*!
introduces a plane `h` that passes through point `p` and 
that has as an orthogonal direction equal to `d`. 
*/ 
Plane_3(const Point_3<Kernel> &p, 
const Direction_3<Kernel>&d); 

/*!
introduces a plane `h` that is defined through the three points 
`l.point(0)`, `l.point(1)` and `p`. 
*/ 
Plane_3(const Line_3<Kernel> &l, 
const Point_3<Kernel> &p); 

/*!
introduces a plane `h` that is defined through the three points 
`r.point(0)`, `r.point(1)` and `p`. 
*/ 
Plane_3(const Ray_3<Kernel> &r, 
const Point_3<Kernel> &p); 

/*!
introduces a plane `h` that is defined through the three points 
`s.source()`, `s.target()` and `p`. 
*/ 
Plane_3(const Segment_3<Kernel> &s, 
const Point_3<Kernel> &p); 

/*!
introduces a plane `h` that is defined as the plane containing 
the circle. 
*/ 
Plane_3(const Circle_3<Kernel> &c); 

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: two planes are equal, iff they have a non 
empty intersection and the same orientation. 
*/ 
bool operator==(const Plane_3<Kernel> &h2) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Plane_3<Kernel> &h2) const; 

/*!
returns the first coefficient of `h`. 
*/ 
Kernel::RT a() const; 

/*!
returns the second coefficient of `h`. 
*/ 
Kernel::RT b() const; 

/*!
returns the third coefficient of `h`. 
*/ 
Kernel::RT c() const; 

/*!
returns the fourth coefficient of `h`. 
*/ 
Kernel::RT d() const; 

/*!
returns the line that is perpendicular to `h` and that 
passes through point `p`. The line is oriented from 
the negative to the positive side of `h`. 
*/ 
Line_3<Kernel> perpendicular_line(const Point_3<Kernel> &p) const; 

/*!
returns the orthogonal projection of `p` on `h`. 
*/ 
Point_3<Kernel> projection(const Point_3<Kernel> &p) const; 

/*!
returns the plane with opposite orientation. 
*/ 
Plane_3<Kernel> opposite() const; 

/*!
returns an arbitrary point on `h`. 
*/ 
Point_3<Kernel> point() const; 

/*!
returns a vector that is orthogonal to `h` and that 
is directed to the positive side of `h`. 
*/ 
Vector_3<Kernel> orthogonal_vector() const; 

/*!
returns the direction that is orthogonal to `h` and that 
is directed to the positive side of `h`. 
*/ 
Direction_3<Kernel> orthogonal_direction() const; 

/*!
returns a vector orthogonal to 
`orthogonal_vector()`. 
*/ 
Vector_3<Kernel> base1() const; 

/*!
returns a vector that is both orthogonal to `base1()`, 
and to `orthogonal_vector()`, and such that the result of 
`orientation( point(), point() + base1(), point()+base2(), point() + orthogonal_vector() )` is positive. 
*/ 
Vector_3<Kernel> base2() const; 

/// @} 

/// \name 2D Conversion 
/// The following functions provide conversion between a plane and
/// \cgal's two-dimensional space. The transformation is affine, but
/// not necessarily an isometry. This means, the transformation
/// preserves combinatorics, but not distances.
/// @{

/*!
returns the image point of the projection of `p` 
under an affine transformation, which maps `h` onto the 
\f$ xy\f$-plane, with the \f$ z\f$-coordinate removed. 
*/ 
Point_2<Kernel> to_2d(const Point_3<Kernel> &p) const; 

/*!
returns a point `q`, such that `to_2d( to_3d( p ))` 
is equal to `p`. 
*/ 
Point_3<Kernel> to_3d(const Point_2<Kernel> &p) const; 

/// @} 

/// \name Predicates 
/// @{

/*!
returns either \ref ON_ORIENTED_BOUNDARY, or 
the constant \ref ON_POSITIVE_SIDE, or the constant 
\ref ON_NEGATIVE_SIDE, 
determined by the position of `p` relative to the oriented plane `h`. 

*/ 
Oriented_side oriented_side(const Point_3<Kernel> &p) const; 

/// @}

/// \name Convenience Boolean Functions
/// @{

/*!

*/ 
bool has_on(const Point_3<Kernel> &p) const; 

/*!

*/ 
bool has_on_positive_side(const Point_3<Kernel> &p) const; 

/*!

*/ 
bool has_on_negative_side(const Point_3<Kernel> &p) const; 

/*!

*/ 
bool has_on(const Line_3<Kernel> &l) const; 

/*!

*/ 
bool has_on(const Circle_3<Kernel> &l) const; 

/*!
Plane `h` is degenerate, if the coefficients `a`, 
`b`, and `c` of the plane equation are zero. 
*/ 
bool is_degenerate() const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
returns the plane obtained by applying `t` on a point of `h` 
and the orthogonal direction of `h`. 
*/ 
Plane_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const; 

/// @}

}; /* end Plane_3 */
} /* end namespace CGAL */
