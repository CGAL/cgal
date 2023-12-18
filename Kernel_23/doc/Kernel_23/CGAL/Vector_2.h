namespace CGAL {

/*!
\ingroup kernel_classes2

An object `v` of the class `Vector_2` is a vector in the two-dimensional
vector space \f$ \mathbb{R}^2\f$. Geometrically spoken, a vector is the difference
of two points \f$ p_2\f$, \f$ p_1\f$ and denotes the direction and the distance
from \f$ p_1\f$ to \f$ p_2\f$.

\cgal defines a symbolic constant \ref NULL_VECTOR. We
will explicitly state where you can pass this constant as an argument
instead of a vector initialized with zeros.

\cgalModels{Kernel::Vector_2,Hashable if `Kernel` is a %Cartesian kernel and if `Kernel::FT` is `Hashable`}

*/
template< typename Kernel >
class Vector_2 {
public:

/// \name Types
/// @{

/*!
An iterator for enumerating the
%Cartesian coordinates of a vector.
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

/// \name Creation
/// @{

/*!
introduces the vector `b-a`.
*/
Vector_2(const Point_2<Kernel> &a, const Point_2<Kernel> &b);

/*!
introduces the vector `s.target()-s.source()`.
*/
Vector_2(const Segment_2<Kernel> &s);

/*!
introduces the vector having the same direction as `r`.
*/
Vector_2(const Ray_2<Kernel> &r);

/*!
introduces the vector having the same direction as `l`.
*/
Vector_2(const Line_2<Kernel> &l);

/*!
introduces a null vector `v`.
*/
Vector_2(const Null_vector &NULL_VECTOR);

/*!
introduces a vector `v` initialized to `(x,y)`.
*/
Vector_2(int x, int y);

/*!
introduces a vector `v` initialized to `(x,y)`.
\cgalEpicExact
*/
Vector_2(double x, double y);

/*!
introduces a vector `v` initialized to `(hx/hw,hy/hw)`.
\pre `hw != 0`.
*/
Vector_2(const Kernel::RT &hx, const Kernel::RT &hy, const Kernel::RT &hw = RT(1));

/*!
introduces a vector `v` initialized to `(x,y)`.
\cgalEpicExact
*/
Vector_2(const Kernel::FT &x, const Kernel::FT &y);

/// @}

/// \name Coordinate Access
/// There are two sets of coordinate access functions, namely to the
/// homogeneous and to the %Cartesian coordinates. They can be used
/// independently from the chosen kernel model. Note that you do not
/// lose information with the homogeneous representation, because the
/// `FieldNumberType` is a quotient.
/// @{

/*!
returns the homogeneous \f$ x\f$ coordinate.
*/
Kernel::RT hx() const;

/*!
returns the homogeneous \f$ y\f$ coordinate.
*/
Kernel::RT hy() const;

/*!
returns the homogenizing coordinate.
*/
Kernel::RT hw() const;

/*!
returns the `x`-coordinate of `v`, that is `hx()`/`hw()`.
\cgalEpicExact
*/
Kernel::FT x() const;

/*!
returns the `y`-coordinate of `v`, that is `hy()`/`hw()`.
\cgalEpicExact
*/
Kernel::FT y() const;

/// @}

/// \name Convenience Operators
/// The following operations are for convenience and for compatibility
/// with higher dimensional vectors. Again they come in a %Cartesian
/// and homogeneous flavor.
/// @{


/*!
returns the i'th homogeneous coordinate of `v`.
\pre `0 <= i <= 2`.
*/
Kernel::RT homogeneous(int i) const;

/*!
returns the i'th %Cartesian coordinate of `v`.
\pre `0 <= i <= 1`.

\cgalEpicExact
*/
Kernel::FT cartesian(int i) const;

/*!
returns `cartesian(i)`.
\pre `0 <= i <= 1`.

\cgalEpicExact
*/
Kernel::FT operator[](int i) const;

/*!
returns an iterator to the %Cartesian coordinates
of `v`, starting with the 0th coordinate.
*/
Cartesian_const_iterator cartesian_begin() const;

/*!
returns an off the end iterator to the %Cartesian
coordinates of `v`.
*/
Cartesian_const_iterator cartesian_end() const;

/*!
returns the dimension (the constant 2).
*/
int dimension() const;

/*!
returns the direction which passes through `v`.
\cgalEpicExact
*/
Direction_2<Kernel> direction() const;

/*!
returns the vector obtained by applying `t` on `v`.
*/
Vector_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const;

/*!
returns the vector perpendicular to `v` in clockwise or
counterclockwise orientation.
*/
Vector_2<Kernel> perpendicular(const Orientation &o) const;

/// @}

/// \name Operators
/// @{

/*!
Test for equality: two vectors are equal, iff their \f$ x\f$ and \f$ y\f$
coordinates are equal. You can compare a vector with the
`NULL_VECTOR`.
*/
bool operator==(const Vector_2<Kernel> &w) const;

/*!
Test for inequality. You can compare a vector with the
`NULL_VECTOR`.
*/
bool operator!=(const Vector_2<Kernel> &w) const;

/*!
Addition.
*/
Vector_2<Kernel> operator+(const Vector_2<Kernel> &w) const;

/*!
Addition.
*/
Vector_2<Kernel>& operator+=(const Vector_2<Kernel> &w);


/*!
Subtraction.
*/
Vector_2<Kernel> operator-(const Vector_2<Kernel> &w) const;

/*!
Subtraction.
*/
Vector_2<Kernel>& operator-=(const Vector_2<Kernel> &w);

/*!
returns the opposite vector.
\cgalEpicExact
*/
Vector_2<Kernel> operator-() const;

/*!
returns the scalar product (= inner product) of the two vectors.
*/
Kernel::FT operator*(const Vector_2<Kernel> &w) const;

/*!
Division by a scalar.
*/
Vector_2<Kernel> operator/(const Kernel::RT &s) const;

/*!
Division by a scalar.
*/
Vector_2<Kernel>& operator/=(const Kernel::RT &s);

/*!
Multiplication by a scalar.
*/
Vector_2<Kernel>& operator*=(const Kernel::RT &s);

/// @}

/*!
returns the squared length of `v`.
*/
Kernel::FT squared_length() const;

}; /* end Vector_2 */


/// \ingroup Kernel_operator_prod
///@{

/*!
Multiplication with a scalar from the right.
\relates Vector_2
*/
Vector_2<Kernel>
operator*(const Vector_2<Kernel> &v, const Kernel::RT &s);

/*!
Multiplication with a scalar from the right.
\relates Vector_2
*/
Vector_2<Kernel>
operator*(const Vector_2<Kernel> &v, const Kernel::FT &s);

/*!
Multiplication with a scalar from the left.
\relates Vector_2
*/
Vector_2<Kernel>
operator*(const Kernel::RT &s, const Vector_2<Kernel> &v);

/*!
Multiplication with a scalar from the left.
\relates Vector_2
*/
Vector_2<Kernel>
operator*(const Kernel::FT &s, const Vector_2<Kernel> &v);

/// @}

} /* end namespace CGAL */
