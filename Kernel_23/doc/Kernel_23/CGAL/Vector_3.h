namespace CGAL {

/*!
\ingroup kernel_classes3

An object of the class `Vector_3` is a vector in the three-dimensional
vector space \f$ \mathbb{R}^3\f$. Geometrically spoken a vector is the difference
of two points \f$ p_2\f$, \f$ p_1\f$ and denotes the direction and the distance
from \f$ p_1\f$ to \f$ p_2\f$.

\cgal defines a symbolic constant \ref NULL_VECTOR. We
will explicitly state where you can pass this constant as an argument
instead of a vector initialized with zeros.

\cgalModels{Kernel::Vector_3,Hashable if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`}

\sa `cross_product_grp`
\sa `determinant_grp`

*/
template< typename Kernel >
class Vector_3 {
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
Vector_3(const Point_3<Kernel> &a, const Point_3<Kernel> &b);

/*!
introduces the vector `s.target()-s.source()`.
*/
Vector_3(const Segment_3<Kernel> &s);

/*!
introduces a vector having the same direction as `r`.
*/
Vector_3(const Ray_3<Kernel> &r);

/*!
introduces a vector having the same direction as `l`.
*/
Vector_3(const Line_3<Kernel> &l);

/*!
introduces a null vector `v`.
\cgalEpicExact
*/
Vector_3(const Null_vector &NULL_VECTOR);

/*!
introduces a vector `v` initialized to `(x, y, z)`.
*/
Vector_3(int x, int y, int z);

/*!
introduces a vector `v` initialized to `(x, y, z)`.
\cgalEpicExact
*/
Vector_3(double x, double y, double z);

/*!
introduces a vector `v` initialized to `(hx/hw, hy/hw, hz/hw)`.
*/
Vector_3(const Kernel::RT &hx, const Kernel::RT &hy, const Kernel::RT &hz, const Kernel::RT &hw = RT(1));

/*!
introduces a vector `v` initialized to `(x, y, z)`.
\cgalEpicExact
*/
Vector_3(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z);

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
returns the homogeneous \f$ z\f$ coordinate.
*/
Kernel::RT hz() const;

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

/*!
returns the `z` coordinate of `v`, that is `hz()`/`hw()`.
\cgalEpicExact
*/
Kernel::FT z() const;

/// @}

/// \name Convenience Operations
/// The following operations are for convenience and for compatibility
/// with higher dimensional vectors. Again they come in a %Cartesian
/// and homogeneous flavor.
/// @{

/*!
returns the i-th homogeneous coordinate of `v`.
\pre `0 <= i <= 3`.
*/
Kernel::RT homogeneous(int i) const;

/*!
returns the i-th %Cartesian coordinate of `v`.
\pre  `0 <= i <= 2`

\cgalEpicExact
*/
Kernel::FT cartesian(int i) const;

/*!
returns `cartesian(i)`.
\pre  `0 <= i <= 2`

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
returns the dimension (the constant 3).
*/
int dimension() const;

/*!
returns the vector obtained by applying `t` on `v`.
*/
Vector_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const;

/*!
returns the direction of `v`.

\cgalEpicExact
*/
Direction_3<Kernel> direction() const;

/// @}

/// \name Operators
/// @{

/*!
Test for equality: two vectors are equal, iff their \f$ x\f$, \f$ y\f$
and \f$ z\f$ coordinates are equal. You can compare a vector with the
`NULL_VECTOR`.
*/
bool operator==(const Vector_3<Kernel> &w) const;

/*!
Test for inequality. You can compare a vector with the
`NULL_VECTOR`.
*/
bool operator!=(const Vector_3<Kernel> &w) const;

/*!
Addition.
*/
Vector_3<Kernel> operator+(const Vector_3<Kernel> &w) const;

/*!
Addition.
*/
Vector_3<Kernel>& operator+=(const Vector_3<Kernel> &w);

/*!
Subtraction.
*/
Vector_3<Kernel> operator-(const Vector_3<Kernel> &w) const;

/*!
Subtraction.
*/
Vector_3<Kernel>& operator-=(const Vector_3<Kernel> &w);

/*!
Returns the opposite vector.
\cgalEpicExact
*/
Vector_3<Kernel> operator-() const;

/*!
Division by a scalar.
*/
Vector_3<Kernel> operator/(const Kernel::RT &s) const;

/*!
Division by a scalar.
*/
Vector_3<Kernel>& operator/=(const Kernel::RT &s);

/*!
returns the scalar product (= inner product) of the two vectors.
*/
Kernel::FT operator*(const Vector_3<Kernel> &w) const;

/*!
Multiplication by a scalar.
*/
Vector_3<Kernel>& operator*=(const Kernel::FT &s);

/// @}

/*!
returns the squared length of `v`.
*/
Kernel::FT squared_length() const;

}; /* end Vector_3 */

/// \ingroup Kernel_operator_prod
///@{

/*!
Multiplication with a scalar from the right.
\relates Vector_3
*/
Vector_3<Kernel>
operator*(const Vector_3<Kernel> &v, const Kernel::RT &s);

/*!
Multiplication with a scalar from the right.
\relates Vector_3
*/
Vector_3<Kernel>
operator*(const Vector_3<Kernel> &v, const Kernel::FT &s);

/*!
Multiplication with a scalar from the left.
\relates Vector_3
*/
Vector_3<Kernel>
operator*(const Kernel::RT &s, const Vector_3<Kernel> &v);

/*!
Multiplication with a scalar from the left.
\relates Vector_3
*/
Vector_3<Kernel>
operator*(const Kernel::FT &s, const Vector_3<Kernel> &v);

/// @}

} /* end namespace CGAL */
