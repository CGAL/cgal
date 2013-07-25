namespace CGAL {
/*!
\ingroup kernel_classes2

An object `l` of the data type `Line_2` is a directed 
straight line in the two-dimensional Euclidean plane \f$ \E^2\f$. It is 
defined by the set of points with %Cartesian coordinates \f$ (x,y)\f$ 
that satisfy the equation 

\f[ l:\; a\, x +b\, y +c = 0. \f]

The line splits \f$ \E^2\f$ in a <I>positive</I> and a <I>negative</I> 
side. A point `p` with %Cartesian coordinates 
\f$ (px, py)\f$ is on the positive side of `l`, iff 
\f$ a\, px + b\, py +c > 0\f$, it is 
on the negative side of `l`, iff 
\f$ a\, px + b\, py +c < 0\f$. 
The positive side is to the left of `l`. 

\cgalHeading{Example}

Let us first define two %Cartesian two-dimensional points in the Euclidean 
plane \f$ \E^2\f$. Their 
dimension and the fact that they are %Cartesian is expressed by 
the suffix `_2` and the representation type `Cartesian`. 

\code
Point_2< Cartesian<double> > p(1.0,1.0), q(4.0,7.0); 
\endcode

To define a line `l` we write: 

\code
Line_2< Cartesian<double> > l(p,q); 
\endcode

\sa `Kernel::Line_2` 

*/
template< typename Kernel >
class Line_2 {
public:

/// \name Creation 
/// @{

/*!
introduces a line `l` with the line equation in %Cartesian 
coordinates \f$ ax +by +c = 0\f$. 
*/ 
Line_2(const Kernel::RT &a, const Kernel::RT &b, const Kernel::RT &c); 

/*!
introduces a line `l` passing through the points `p` and `q`. 
Line `l` is directed from `p` to `q`. 
*/ 
Line_2(const Point_2<Kernel> &p, const Point_2<Kernel> &q); 

/*!
introduces a line `l` passing through point `p` with 
direction `d`. 
*/ 
Line_2(const Point_2<Kernel> &p, const Direction_2<Kernel>&d); 

/*!
introduces a line `l` passing through point `p` and 
oriented by `v`. 
*/ 
Line_2(const Point_2<Kernel> &p, const Vector_2<Kernel>&v); 

/*!
introduces a line `l` supporting the segment `s`, 
oriented from source to target. 
*/ 
Line_2(const Segment_2<Kernel> &s); 

/*!
introduces a line `l` supporting the ray `r`, 
with same orientation. 
*/ 
Line_2(const Ray_2<Kernel> &r); 

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: two lines are equal, iff they have a non 
empty intersection and the same direction. 
*/ 
bool operator==(const Line_2<Kernel> &h) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Line_2<Kernel> &h) const; 

/*!
returns the first coefficient of `l`. 
*/ 
Kernel::RT a() const; 

/*!
returns the second coefficient of `l`. 
*/ 
Kernel::RT b() const; 

/*!
returns the third coefficient of `l`. 
*/ 
Kernel::RT c() const; 

/*!
returns an arbitrary point on `l`. It holds 
`point(i) == point(j)`, iff `i==j`. 
Furthermore, `l` is directed from `point(i)` 
to `point(j)`, for all `i` \f$ <\f$ `j`. 
*/ 
Point_2<Kernel> point(int i) const; 

/*!
returns the orthogonal projection of `p` onto `l`. 
*/ 
Point_2<Kernel> projection(const Point_2<Kernel> &p) const; 

/*!
returns the \f$ x\f$-coordinate of the point at `l` with 
given \f$ y\f$-coordinate. 
\pre `l` is not horizontal. 
*/ 
Kernel::FT x_at_y(const Kernel::FT &y) const; 

/*!
returns the \f$ y\f$-coordinate of the point at `l` with 
given \f$ x\f$-coordinate. 
\pre `l` is not vertical. 
*/ 
Kernel::FT y_at_x(const Kernel::FT &x) const; 

/// @} 

/// \name Predicates 
/// @{

/*!
line `l` is degenerate, if the coefficients `a` and 
`b` of the line equation are zero. 
*/ 
bool is_degenerate() const; 

/*!

*/ 
bool is_horizontal() const; 

/*!

*/ 
bool is_vertical() const; 

/*!
returns \ref ON_ORIENTED_BOUNDARY, 
\ref ON_NEGATIVE_SIDE, or the constant 
\ref ON_POSITIVE_SIDE, 
depending on the position of `p` relative to the oriented line `l`. 

*/ 
Oriented_side oriented_side(const Point_2<Kernel> &p) const; 

/// @}

/// \name Convenience Boolean Functions
/// @{

/*!

*/ 
bool has_on(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_positive_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_negative_side(const Point_2<Kernel> &p) const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
returns a vector having the direction of `l`. 
*/ 
Vector_2<Kernel> to_vector() const; 

/*!
returns the direction of `l`. 
*/ 
Direction_2<Kernel> direction() const; 

/*!
returns the line with opposite direction. 
*/ 
Line_2<Kernel> opposite() const; 

/*!
returns the line perpendicular to `l` and passing through `p`, 
where the direction is the direction of `l` rotated 
counterclockwise by 90 degrees. 
*/ 
Line_2<Kernel> perpendicular(const Point_2<Kernel> &p) const; 

/*!
returns the line obtained by applying `t` on a point on `l` 
and the direction of `l`. 
*/ 
Line_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const; 

/// @}

}; /* end Line_2 */
} /* end namespace CGAL */

