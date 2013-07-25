namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance \f$ s\f$ of the data type `Segment_d` is a directed 
straight line segment in \f$ d\f$-dimensional Euclidean space connecting 
two points \f$ p\f$ and \f$ q\f$. \f$ p\f$ is called the source point and \f$ q\f$ is 
called the target point of \f$ s\f$, both points are called endpoints of 
\f$ s\f$. A segment whose endpoints are equal is called <I>degenerate</I>. 

\cgalHeading{Implementation}

Segments are implemented by a pair of points as an item type. All 
operations like creation, initialization, tests, the calculation of 
the direction and source - target vector, input and output on a 
segment \f$ s\f$ take time \f$ O(s.dimension())\f$. `dimension()`, 
coordinate and end point access, and identity test take constant time. 
The operations for intersection calculation also take time 
\f$ O(s.dimension())\f$. The space requirement is 
\f$ O(s.dimension())\f$. 

*/
template< typename Kernel >
class Segment_d {
public:

/// \name Types 
/// @{

/*!
the linear algebra layer. 
*/ 
typedef unspecified_type LA; 

/// @} 

/// \name Creation 
/// @{

/*!
introduces a variable `s` 
of type `Segment_d<Kernel>`. 
*/ 
Segment_d<Kernel>(); 

/*!
introduces a 
variable `s` of type `Segment_d<Kernel>` which is initialized to 
the segment \f$ (p,q)\f$.

\pre `p.dimension()==q.dimension()`. 
*/ 
Segment_d<Kernel>(Point_d<Kernel> p, Point_d<Kernel> q); 

/*!
introduces a 
variable `s` of type `Segment_d<Kernel>` which is initialized to 
the segment `(p,p+v)`.

\pre `p.dimension()==v.dimension()`. 
*/ 
Segment_d<Kernel>(Point_d<Kernel> p, Vector_d<Kernel> v); 

/// @} 

/// \name Operations 
/// @{

/*!
returns the dimension of the ambient 
space. 
*/ 
int dimension(); 

/*!
returns the source point of segment 
`s`. 
*/ 
Point_d<Kernel> source(); 

/*!
returns the target point of segment 
`s`. 
*/ 
Point_d<Kernel> target(); 

/*!
returns source or target of 
`s`: `vertex(0)` returns the source, `vertex(1)` returns 
the target. The parameter \f$ i\f$ is taken modulo \f$ 2\f$, which gives easy 
access to the other vertex.

\pre \f$ i \geq0\f$. 
*/ 
Point_d<Kernel> vertex(int i) ; 

/*!
returns `vertex(i)`. 
*/ 
Point_d<Kernel> point(int i) ; 

/*!
returns `vertex(i)`. 
*/ 
Point_d<Kernel> operator[](int i) ; 

/*!
returns the lexicographically smaller 
vertex. 
*/ 
Point_d<Kernel> min() ; 

/*!
returns the lexicographically larger 
vertex. 
*/ 
Point_d<Kernel> max() ; 

/*!
returns the segment 
`(target(),source())`. 
*/ 
Segment_d<Kernel> opposite() ; 

/*!
returns the direction from 
source to target.

\pre `s` is non-degenerate. 
*/ 
Direction_d<Kernel> direction() ; 

/*!
returns the vector from source to 
target. 
*/ 
Vector_d<Kernel> vector() ; 

/*!
returns the square of the length of 
`s`. 
*/ 
FT squared_length() ; 

/*!
returns true if \f$ p\f$ lies 
on `s` and false otherwise.

\pre `s.dimension()==p.dimension()`. 
*/ 
bool has_on(const Point_d<Kernel>& p) ; 

/*!
returns the supporting line 
of `s`.

\pre `s` is non-degenerate. 
*/ 
Line_d<Kernel> supporting_line() ; 

/*!
returns \f$ t(s)\f$.

\pre `s.dimension()==t.dimension()`. 
*/ 
Segment_d<Kernel> transform(const Aff_transformation_d<Kernel>& t) 
; 

/*!
returns 
\f$ s+v\f$, i.e., `s` translated by vector \f$ v\f$.

\pre `s.dimension()==v.dimension()`. 
*/ 
Segment_d<Kernel> operator+(const Vector_d<Kernel>& v) ; 

/*!
returns true if `s` is 
degenerate i.e.\ `s.source()=s.target()`. 
*/ 
bool is_degenerate() ; 

/// @}

}; /* end Segment_d */

/*!
Test for equality as unoriented 
segments.

\pre `s1.dimension()==s2.dimension()`. 
\relates Segment_d 
*/ 
bool weak_equality(const Segment_d<Kernel>& s1, const Segment_d<Kernel>& s2) ; 

/*!
return true if one of the segments is degenerate or if the 
unoriented supporting lines are parallel.

\pre `s1.dimension()==s2.dimension()`. 
\relates Segment_d 
*/ 
bool parallel(const Segment_d<Kernel>& s1, const Segment_d<Kernel>& s2) ; 


/*!
if `s1` and `s2` 
touch in a common end point, this point is assigned to `common` 
and the result is `true`, otherwise the result is `false`. If 
`s1==s2` then one of the endpoints is returned.

\pre `s1.dimension()==s2.dimension()`. 
\relates Segment_d 
*/ 
bool common_endpoint(const Segment_d<Kernel>& s1, const Segment_d<Kernel>& s2, Point_d<Kernel>& common) ; 

} /* end namespace CGAL */
