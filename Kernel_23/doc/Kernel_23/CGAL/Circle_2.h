namespace CGAL {

/*!
\ingroup kernel_classes2

An object `c` of type `Circle_2` is a circle in the 
two-dimensional Euclidean plane \f$ \E^2\f$. The circle is oriented, i.e.\ its 
boundary has clockwise or counterclockwise orientation. The 
boundary splits \f$ \E^2\f$ into a positive and a negative side, where the 
positive side is to the left of the boundary. The boundary also 
splits \f$ \E^2\f$ into a bounded and an unbounded side. Note that the 
circle can be degenerated, i.e.\ the squared radius may be zero. 

\sa `Kernel::Circle_2` 

*/
template< typename Kernel >
class Circle_2 {
public:

/// \name Creation 
/// @{

/*!

introduces a variable `c` of type `Circle_2`. 
It is initialized to the circle with center `center`, 
squared radius `squared_radius` and orientation 
`ori`. 
\pre `ori` \f$ \neq\f$ `COLLINEAR`, and further, `squared_radius` \f$ \geq\f$ 0. 
*/ 
Circle_2(const Point_2<Kernel> &center, 
const Kernel::FT &squared_radius, 
const Orientation &ori = COUNTERCLOCKWISE); 

/*!

introduces a variable `c` of type `Circle_2`. 
It is initialized to the unique circle which passes through 
the points `p`, `q` and `r`. The orientation of 
the circle is the orientation of the point triple `p`, 
`q`, `r`. 
\pre `p`, `q`, and `r` are not collinear. 
*/ 
Circle_2(const Point_2<Kernel> &p, 
const Point_2<Kernel> &q, 
const Point_2<Kernel> &r); 

/*!

introduces a variable `c` of type `Circle_2`. 
It is initialized to the circle with diameter \f$ \overline{pq}\f$ 
and orientation `ori`. 
\pre `ori` \f$ \neq\f$ `COLLINEAR`. 
*/ 
Circle_2( const Point_2<Kernel> &p, 
const Point_2<Kernel> &q, 
const Orientation &ori = COUNTERCLOCKWISE); 

/*!

introduces a variable `c` of type `Circle_2`. 
It is initialized to the circle with center `center`, squared 
radius zero and orientation `ori`. 
\pre `ori` \f$ \neq\f$ `COLLINEAR`. 
\post `c.is_degenerate()` = `true`. 
*/ 
Circle_2( const Point_2<Kernel> &center, 
          const Orientation &ori = COUNTERCLOCKWISE); 

/// @} 

/// \name Access Functions 
/// @{

/*!

returns the center of `c`. 
*/ 
const Point_2<Kernel> &center( ) const; 

/*!

returns the squared radius of `c`. 
*/ 
const Kernel::FT& squared_radius( ) const; 

/*!

returns the orientation of `c`. 
*/ 
Orientation orientation( ) const; 

/*!

returns `true`, iff `c` and `circle2` are equal, 
i.e.\ if they have the same center, same squared radius and 
same orientation. 
*/ 
bool operator == ( const Circle_2<Kernel>& circle2) const; 

/*!

returns `true`, iff `c` and `circle2` are not equal. 
*/ 
bool operator != ( const Circle_2<Kernel> & circle2) const; 

/// @} 

/// \name Predicates 
/// @{

/*!

returns `true`, iff `c` is degenerate, i.e.\ if `c` has squared radius zero. 
*/ 
bool is_degenerate( ) const; 

/*!

returns either the constant \ref ON_ORIENTED_BOUNDARY, 
\ref ON_POSITIVE_SIDE, or \ref ON_NEGATIVE_SIDE,
iff `p` lies on the boundary, properly on the 
positive side, or properly on the negative side 
of `c`, resp. 
*/ 
Oriented_side 
oriented_side( const Point_2<Kernel> &p) const; 

/*!

returns \ref ON_BOUNDED_SIDE, 
\ref ON_BOUNDARY, or \ref ON_UNBOUNDED_SIDE 
iff `p` lies properly inside, on the boundary, or properly 
outside of `c`, resp. 
*/ 
Bounded_side 
bounded_side( const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_positive_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_negative_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_boundary(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_bounded_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_unbounded_side(const Point_2<Kernel> &p) const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!

returns the circle with the same center and squared radius as 
`c` but with opposite orientation. 
*/ 
Circle_2<Kernel> opposite() const; 

/*!

returns the circle obtained by applying \f$ at\f$ on `c`. 
\pre `at` is an orthogonal transformation. 
*/ 
Circle_2<Kernel> orthogonal_transform( 
Aff_transformation_2<Kernel> const& at) const; 

/*!

returns a bounding box containing `c`. 
*/ 
Bbox_2 bbox() const; 

/// @}

}; /* end Circle_2 */
} /* end namespace CGAL */
