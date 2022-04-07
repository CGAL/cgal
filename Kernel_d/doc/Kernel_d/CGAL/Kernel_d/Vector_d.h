namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of data type `Vector_d<Kernel>` is a vector of Euclidean
space in dimension \f$ d\f$. A vector \f$ r = (r_0,\ldots,r_{ d - 1})\f$ can be
represented in homogeneous coordinates \f$ (h_0,\ldots,h_d)\f$ of number
type `RT`, such that \f$ r_i = h_i/h_d\f$ which is of type `FT`. We
call the \f$ r_i\f$'s the %Cartesian coordinates of the vector. The
homogenizing coordinate \f$ h_d\f$ is positive.

This data type is meant for use in computational geometry. It realizes
free vectors as opposed to position vectors (type `Point_d`). The
main difference between position vectors and free vectors is their
behavior under affine transformations, e.g., free vectors are
invariant under translations.

\cgalHeading{Downward compatibility}

We provide all operations of the
lower dimensional interface `x()`, `y()`, `z()`,
`hx()`, `hy()`, `hz()`, `hw()`.

\cgalHeading{Implementation}

Vectors are implemented by arrays of variables of type `RT`. All
operations like creation, initialization, tests, vector arithmetic,
input and output on a vector \f$ v\f$ take time \f$ O(v.dimension())\f$.
coordinate access, `dimension()` and conversions take constant
time. The space requirement of a vector is \f$ O(v.dimension())\f$.

*/
template< typename Kernel >
class Vector_d {
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

/*!
construction tag.
*/
typedef unspecified_type Base_vector;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable `v`
of type `Vector_d<Kernel>`.
*/
Vector_d<Kernel>();

/*!
introduces the zero
vector `v` of type `Vector_d<Kernel>` in \f$ d\f$-dimensional space.
For the creation flag `CGAL::NULL_VECTOR` can be used.
*/
Vector_d<Kernel>(int d, Null_vector);

/*!
introduces a variable
`v` of type `Vector_d<Kernel>` in dimension `d`. If
`size [first,last) == d` this creates a vector with %Cartesian
coordinates `set [first,last)`. If `size [first,last) == p+1` the range specifies the homogeneous coordinates \f$ H = set
[first,last) = (\pm h_0, \pm h_1, \ldots, \pm h_d)\f$ where the
sign chosen is the sign of \f$ h_d\f$.

\pre `d` is nonnegative, `[first,last)` has `d` or `d+1` elements where the last has to be non-zero.
\tparam InputIterator has `RT` as value type.
*/
template <class InputIterator> Vector_d<Kernel>(int d,
InputIterator first, InputIterator last);

/*!
introduces a
variable `v` of type `Vector_d<Kernel>` in dimension `d`
initialized to the vector with homogeneous coordinates as defined by
`H = set [first,last)` and `D`: \f$ (\pm H_0,
\pm H_1, \ldots, \pm H_{D-1}, \pm H_D)\f$. The sign
chosen is the sign of \f$ D\f$.

\pre `D` is non-zero, the iterator range defines a \f$ d\f$-tuple of `RT`.
\tparam InputIterator has `RT` as value type.
*/
template <class InputIterator> Vector_d<Kernel>(int d,
InputIterator first, InputIterator last, RT D);

/*!
returns a
variable `v` of type `Vector_d<Kernel>` initialized to the \f$ i\f$-th
base vector of dimension \f$ d\f$.

\pre \f$ 0 \leq i < d\f$.
*/
Vector_d<Kernel>(int d, Base_vector, int i);

/*!
introduces a
variable `v` of type `Vector_d<Kernel>` in \f$ 2\f$-dimensional space.
\pre \f$ w \neq0\f$.
*/
Vector_d<Kernel>(RT x, RT y, RT w = 1);

/*!
introduces a
variable `v` of type `Vector_d<Kernel>` in \f$ 3\f$-dimensional space.
\pre \f$ w \neq0\f$.
*/
Vector_d<Kernel>(RT x, RT y, RT z, RT w);

/// @}

/// \name Operations
/// @{

/*!
returns the dimension of `v`.
*/
int dimension() ;

/*!
returns the \f$ i\f$-th %Cartesian
coordinate of `v`.

\pre \f$ 0 \leq i < d\f$.
*/
FT cartesian(int i) ;

/*!
returns the \f$ i\f$-th %Cartesian
coordinate of `v`.

\pre \f$ 0 \leq i < d\f$.
*/
FT operator[](int i) ;

/*!
returns the \f$ i\f$-th homogeneous
coordinate of `v`.

\pre \f$ 0 \leq i \leq d\f$.
*/
RT homogeneous(int i) ;

/*!
returns the square of the length of
`v`.
*/
FT squared_length() ;

/*!
returns an
iterator pointing to the zeroth %Cartesian coordinate of `v`.
*/
Cartesian_const_iterator cartesian_begin() ;

/*!
returns an
iterator pointing beyond the last %Cartesian coordinate of `v`.

*/
Cartesian_const_iterator cartesian_end() ;

/*!
returns an
iterator pointing to the zeroth homogeneous coordinate of `v`.

*/
Homogeneous_const_iterator homogeneous_begin() ;

/*!
returns an
iterator pointing beyond the last homogeneous coordinate of `v`.

*/
Homogeneous_const_iterator homogeneous_end() ;

/*!
returns the direction of
`v`.
*/
Direction_d<Kernel> direction() ;

/*!
returns \f$ t(v)\f$.
*/
Vector_d<Kernel> transform(const Aff_transformation_d<Kernel>& t)
;

/// @}

/// \name Arithmetic Operators, Tests and IO
/// @{

/*!
multiplies all
%Cartesian coordinates by `n`.
*/
Vector_d<Kernel>& operator*=(const RT& n) ;

/*!
multiplies all
%Cartesian coordinates by `r`.
*/
Vector_d<Kernel>& operator*=(const FT& r) ;

/*!
returns the vector
with %Cartesian coordinates \f$ v_i/n, 0 \leq i < d\f$.
*/
Vector_d<Kernel> operator/(const RT& n) ;

/*!
returns the vector
with %Cartesian coordinates \f$ v_i/r, 0 \leq i < d\f$.
*/
Vector_d<Kernel> operator/(const FT& r) ;

/*!
divides all
%Cartesian coordinates by `n`.
*/
Vector_d<Kernel>& operator/=(const RT& n) ;

/*!
divides all
%Cartesian coordinates by `r`.
*/
Vector_d<Kernel>& operator/=(const FT& r) ;

/*!
inner product, i.e.,
\f$ \sum_{ 0 \le i < d } v_i w_i\f$, where \f$ v_i\f$ and \f$ w_i\f$ are the
%Cartesian coordinates of \f$ v\f$ and \f$ w\f$ respectively.
*/
FT operator* (const Vector_d<Kernel>& w) ;

/*!
returns the
vector with %Cartesian coordinates \f$ v_i+w_i, 0 \leq i < d\f$.
*/
Vector_d<Kernel> operator+(const Vector_d<Kernel>& w) ;

/*!
addition
plus assignment.
*/
Vector_d<Kernel>& operator+=(const Vector_d<Kernel>& w) ;

/*!
returns the
vector with %Cartesian coordinates \f$ v_i-w_i, 0 \leq i < d\f$.
*/
Vector_d<Kernel> operator-(const Vector_d<Kernel>& w) ;

/*!
subtraction
plus assignment.
*/
Vector_d<Kernel>& operator-=(const Vector_d<Kernel>& w) ;

/*!
returns the vector in opposite
direction.
*/
Vector_d<Kernel> operator-() ;

/*!
returns true if `v` is the zero
vector.
*/
bool is_zero() ;

/*!
returns the vector with %Cartesian coordinates \f$ n v_i\f$.
*/
Vector_d<Kernel> operator*(const RT& n, const Vector_d<Kernel>& v);

/*!
returns the vector with %Cartesian coordinates \f$ r v_i, 0 \leq i <
d\f$.
*/
Vector_d<Kernel> operator*(const FT& r, const Vector_d<Kernel>& v);

/// @}

}; /* end Vector_d */
} /* end namespace CGAL */
