namespace CGAL {

/*!
\ingroup kernel_classes3

An object `l` of the data type `Line_3` is a directed 
straight line in the three-dimensional Euclidean space \f$ \E^3\f$. 

\sa `Kernel::Line_3` 

*/
template< typename Kernel >
class Line_3 {
public:

/// \name Creation 
/// @{

/*!
introduces a line `l` passing through the points `p` and `q`. 
Line `l` is directed from `p` to `q`. 
*/ 
Line_3(const Point_3<Kernel> &p, const Point_3<Kernel> &q); 

/*!
introduces a line `l` passing through point `p` with 
direction `d`. 
*/ 
Line_3(const Point_3<Kernel> &p, const Direction_3<Kernel>&d); 

/*!
introduces a line `l` passing through point `p` and 
oriented by `v`. 
*/ 
Line_3(const Point_3<Kernel> &p, const Vector_3<Kernel>&v); 

/*!
returns the line supporting the segment `s`, 
oriented from source to target. 
*/ 
Line_3(const Segment_3<Kernel> &s); 

/*!
returns the line supporting the ray `r`, with the 
same orientation. 
*/ 
Line_3(const Ray_3<Kernel> &r); 

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: two lines are equal, iff they have a non 
empty intersection and the same direction. 
*/ 
bool operator==(const Line_3<Kernel> &h) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Line_3<Kernel> &h) const; 

/*!
returns the orthogonal projection of `p` on `l`. 
*/ 
Point_3<Kernel> projection(const Point_3<Kernel> &p) const; 

/*!
returns an arbitrary point on `l`. It holds 
`point(i) = point(j)`, iff `i=j`. 
*/ 
Point_3<Kernel> point(int i) const; 

/// @} 

/// \name Predicates 
/// @{

/*!
returns `true` iff line `l` is degenerated to a point. 
*/ 
bool is_degenerate() const; 

/*!
returns `true` iff `p` lies on `l`. 
*/ 
bool has_on(const Point_3<Kernel> &p) const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
returns the plane perpendicular to `l` passing through `p`. 
*/ 
Plane_3<Kernel> perpendicular_plane(const Point_3<Kernel> &p) const; 

/*!
returns the line with opposite direction. 
*/ 
Line_3<Kernel> opposite() const; 

/*!
returns a vector having the same direction as `l`. 
*/ 
Vector_3<Kernel> to_vector() const; 

/*!
returns the direction of `l`. 
*/ 
Direction_3<Kernel> direction() const; 

/*!
returns the line obtained by applying `t` on a point on `l` 
and the direction of `l`. 
*/ 
Line_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const; 

/// @}

}; /* end Line_3 */
} /* end namespace CGAL */
