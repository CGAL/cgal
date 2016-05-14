namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of data type `Point_d<Kernel>` is a point of Euclidean 
space in dimension \f$ d\f$. A point \f$ p = (p_0,\ldots,p_{ d - 1 })\f$ in 
\f$ d\f$-dimensional space can be represented by homogeneous coordinates 
\f$ (h_0,h_1,\ldots,h_d)\f$ of number type `RT` such that \f$ p_i = 
h_i/h_d\f$, which is of type `FT`. The homogenizing coordinate \f$ h_d\f$ 
is positive. 

We call \f$ p_i\f$, \f$ 0 \leq i < d\f$ the \f$ i\f$-th %Cartesian coordinate and 
\f$ h_i\f$, \f$ 0 \le i \le d\f$, the \f$ i\f$-th homogeneous coordinate. We call \f$ d\f$ 
the dimension of the point. 

\cgalHeading{Downward compatibility}

We provide operations of the lower 
dimensional interface `x()`, `y()`, `z()`, `hx()`, 
`hy()`, `hz()`, `hw()`. 

\cgalHeading{Implementation}

Points are implemented by arrays of `RT` items. All operations 
like creation, initialization, tests, point - vector arithmetic, input 
and output on a point \f$ p\f$ take time \f$ O(p.dimension())\f$. 
`dimension()`, coordinate access and conversions take constant 
time. The space requirement for points is \f$ O(p.dimension())\f$. 

*/
template< typename Kernel >
class Point_d {
public:

/// \name Types 
/// @{

/*!
the linear algebra layer. 
*/ 
typedef unspecified_type LA; 

/*!
a read-only iterator for the 
%Cartesian coordinates. 
*/ 
typedef unspecified_type Cartesian_const_iterator; 

/*!
a read-only iterator for the 
homogeneous coordinates. 
*/ 
typedef unspecified_type Homogeneous_const_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
introduces a variable `p` of 
type `Point_d<Kernel>`. 
*/ 
Point_d<Kernel>(); 

/*!
introduces a variable 
`p` of type `Point_d<Kernel>` in \f$ d\f$-dimensional space, 
initialized to the origin. 
*/ 
Point_d<Kernel>(int d, Origin); 

/*!
introduces a variable 
`p` of type `Point_d<Kernel>` in dimension `d`. If `size [first,last) == d` this creates a point with %Cartesian coordinates 
`set [first,last)`. If `size [first,last) == d+1` the range 
specifies the homogeneous coordinates 
\f$ H = set [first,last) = (\pm h_0, \pm h_1, \ldots, \pm h_d)\f$
where the sign chosen is the sign of \f$ h_d\f$.

\pre `d` is nonnegative, `[first,last)` has `d` or `d+1` elements where the last has to be non-zero. 
\tparam InputIterator has `RT` as value type.
*/ 
template <class InputIterator> Point_d<Kernel>(int d, 
InputIterator first, InputIterator last); 

/*!
introduces a 
variable `p` of type `Point_d<Kernel>` in dimension `d` 
initialized to the point with homogeneous coordinates as defined by 
`H = set [first,last)` and `D`: \f$ (\pm H[0], 
\pm H[1], \ldots, \pm H[d-1], \pm D)\f$. 
The sign  chosen is the sign of \f$ D\f$.

\pre `D` is non-zero, the iterator range defines a \f$ d\f$-tuple of `RT`. 
\tparam InputIterator has `RT` as value type.
*/ 
template <class InputIterator> Point_d<Kernel>(int d, 
InputIterator first, InputIterator last, RT D); 

/*!
introduces a variable 
`p` of type `Point_d<Kernel>` in \f$ 2\f$-dimensional space. 
\pre \f$ w \neq0\f$. 
*/ 
Point_d<Kernel>(RT x, RT y, RT w = 1); 

/*!
introduces a 
variable `p` of type `Point_d<Kernel>` in \f$ 3\f$-dimensional 
space.

\pre \f$ w \neq0\f$. 
*/ 
Point_d<Kernel>(RT x, RT y, RT z, RT w); 

/// @} 

/// \name Operations 
/// @{

/*!
returns the dimension of `p`. 

*/ 
int dimension() ; 

/*!
returns the \f$ i\f$-th %Cartesian 
coordinate of `p`.

\pre \f$ 0 \leq i < d\f$. 
*/ 
FT cartesian(int i) ; 

/*!
returns the \f$ i\f$-th %Cartesian 
coordinate of `p`.

\pre \f$ 0 \leq i < d\f$. 
*/ 
FT operator[](int i) ; 

/*!
returns the \f$ i\f$-th homogeneous 
coordinate of `p`.

\pre \f$ 0 \leq i \leq d\f$. 
*/ 
RT homogeneous(int i) ; 

/*!
returns an 
iterator pointing to the zeroth %Cartesian coordinate \f$ p_0\f$ of 
`p`. 
*/ 
Cartesian_const_iterator cartesian_begin() ; 

/*!
returns an 
iterator pointing beyond the last %Cartesian coordinate of `p`. 

*/ 
Cartesian_const_iterator cartesian_end() ; 

/*!
returns 
an iterator pointing to the zeroth homogeneous coordinate \f$ h_0\f$ of 
`p`. 
*/ 
Homogeneous_const_iterator homogeneous_begin() ; 

/*!
returns an 
iterator pointing beyond the last homogeneous coordinate of 
`p`. 
*/ 
Homogeneous_const_iterator homogeneous_end() ; 

/*!
returns \f$ t(p)\f$. 
*/ 
Point_d<Kernel> transform(const Aff_transformation_d<Kernel>& t) 
; 

/// @} 

/// \name Arithmetic Operators, Tests and IO 
/// @{

/*!
returns the 
vector \f$ p-O\f$. 
*/ 
Vector_d<Kernel> operator-(const Origin& o) ; 

/*!
returns \f$ p - 
q\f$.

\pre `p.dimension() == q.dimension()`. 
*/ 
Vector_d<Kernel> operator-(const Point_d<Kernel>& q) ; 

/*!
returns \f$ p + 
v\f$.

\pre `p.dimension() == v.dimension()`. 
*/ 
Point_d<Kernel> operator+(const Vector_d<Kernel>& v) ; 

/*!
returns \f$ p - 
v\f$.

\pre `p.dimension() == v.dimension()`. 
*/ 
Point_d<Kernel> operator-(const Vector_d<Kernel>& v) ; 

/*!
adds `v` 
to `p`.

\pre `p.dimension() == v.dimension()`. 
*/ 
Point_d<Kernel>& operator+=(const Vector_d<Kernel>& v); 

/*!
subtracts 
`v` from `p`.

\pre `p.dimension() == v.dimension()`. 
*/ 
Point_d<Kernel>& operator-=(const Vector_d<Kernel>& v); 

/*!
returns true if `p` is 
the origin. 
*/ 
bool operator==(const Origin&) ; 

/*!
returns true iff `p` is lexicographically
smaller than `q` with respect to %Cartesian
lexicographic order of points.
\pre `p.dimension() == q.dimension()`.
*/
bool operator<(const Point_d<Kernel>& q);

/*!
returns true iff `p` is lexicographically
greater than `q` with respect to %Cartesian
lexicographic order of points.
\pre `p.dimension() == q.dimension()`.
*/
bool operator>(const Point_d<Kernel>& q);

/*!
returns true iff `p` is lexicographically
smaller than or equal to `q` with respect to %Cartesian
lexicographic order of points.
\pre `p.dimension() == q.dimension()`.
*/
bool operator<=(const Point_d<Kernel>& q);

/*!
returns true iff `p` is lexicographically
greater than or equal to `q` with respect to %Cartesian
lexicographic order of points.
\pre `p.dimension() == q.dimension()`.
*/
bool operator>=(const Point_d<Kernel>& q);

/// @}

}; /* end Point_d */
} /* end namespace CGAL */
