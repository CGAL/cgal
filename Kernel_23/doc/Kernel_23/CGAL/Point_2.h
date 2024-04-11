namespace CGAL {

/*!
\ingroup kernel_classes2

An object `p` of the class `Point_2` is a point in the two-dimensional
Euclidean plane \f$ \E^2\f$.

Remember that `Kernel::RT` and `Kernel::FT` denote a
`RingNumberType` and a `FieldNumberType`, respectively. For the kernel
model `Cartesian<NT>`, the two types are the same. For the
kernel model `Homogeneous<NT>`, `Kernel::RT` is equal
to `NT`, and `Kernel::FT` is equal to `Quotient<NT>`.

\cgalHeading{Example}

The following declaration creates two points with
%Cartesian double coordinates.

\code
Point_2< Cartesian<double> > p, q(1.0, 2.0);
\endcode

The variable `p` is uninitialized and should first be used on
the left hand side of an assignment.

\code
p = q;

std::cout << p.x() << " " << p.y() << std::endl;
\endcode

\cgalModels{Kernel::Point_2,Hashable if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`}

*/
template< typename Kernel >
class Point_2 {
public:

/// \name Types
/// @{

/*!
An iterator for enumerating the
%Cartesian coordinates of a point.
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable `p` with %Cartesian coordinates
\f$ (0,0)\f$.
\cgalEpicExact
*/
Point_2(const Origin &ORIGIN);

/*!
introduces a point `p` initialized to `(x,y)`.
*/
Point_2(int x, int y);

/*!
introduces a point `p` initialized to `(x,y)`
provided `RT` supports construction from `double`.
\cgalEpicExact
*/
Point_2(double x, double y);

/*!
introduces a point `p` initialized to `(hx/hw,hy/hw)`.
\pre `hw != Kernel::RT(0)`.
*/
Point_2(const Kernel::RT &hx, const Kernel::RT &hy, const Kernel::RT &hw = RT(1));

/*!
introduces a point `p` initialized to `(x,y)`.
\cgalEpicExact
*/
Point_2(const Kernel::FT &x, const Kernel::FT &y);

/*!
introduces a point from a weighted point.
\cgalEpicExact

\warning The `explicit` keyword is used to avoid accidental implicit conversions
         between Point_2 and Weighted_point_2.
*/
explicit Point_2(const Kernel::Weighted_point_2 &wp);

/// @}

/// \name Operations
/// @{

/*!
Test for equality. Two points are equal, iff their \f$ x\f$ and \f$ y\f$
coordinates are equal. The point can be compared with `ORIGIN`.
*/
bool operator==(const Point_2<Kernel> &q) const;

/*!
Test for inequality. The point can be compared with `ORIGIN`.
*/
bool operator!=(const Point_2<Kernel> &q) const;

/*!
translates the point by the vector `v`.
*/
Point_2<Kernel>& operator+=(const Vector_2<Kernel> &v);

/*!
translates the point by the vector -`v`.
*/
Point_2<Kernel>& operator-=(const Vector_2<Kernel> &v);

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
returns the %Cartesian \f$ x\f$ coordinate, that is `hx()`/`hw()`.
\cgalEpicExact
*/
Kernel::FT x() const;

/*!
returns the %Cartesian \f$ y\f$ coordinate, that is `hy()`/`hw()`.
\cgalEpicExact
*/
Kernel::FT y() const;

/// @}

/// \name Convenience Operations
/// The following operations are for convenience and for compatibility
/// with higher dimensional points. Again they come in a %Cartesian and
/// in a homogeneous flavor.
/// @{

/*!
returns the i'th homogeneous coordinate of `p`.
\pre `0 <= i <= 2`.
*/
Kernel::RT homogeneous(int i) const;

/*!
returns the i'th %Cartesian coordinate of `p`.
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
of `p`, starting with the 0th coordinate.
*/
Cartesian_const_iterator cartesian_begin() const;

/*!
returns an off the end iterator to the Cartesian
coordinates of `p`.
*/
Cartesian_const_iterator cartesian_end() const;

/*!
returns the dimension (the constant 2).
*/
int dimension() const;

/*!
returns a bounding box containing `p`.
\cgalEpicExact
*/
Bbox_2 bbox() const;

/*!
returns the point obtained by applying `t` on `p`.
*/
Point_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const;

/// @}

}; /* end Point_2 */


/*!
returns true iff `p` is lexicographically smaller than `q`,
i.e.\ either if `p.x() < q.x()` or if `p.x() == q.x()` and
`p.y() < q.y()`.
\relates Point_2
*/
bool operator<(const Point_2<Kernel> &p,
const Point_2<Kernel> &q);

/*!
returns true iff `p` is lexicographically greater than `q`.
\relates Point_2
*/
bool operator>(const Point_2<Kernel> &p,
const Point_2<Kernel> &q);

/*!
returns true iff `p` is lexicographically smaller or equal to `q`.
\relates Point_2
*/
bool operator<=(const Point_2<Kernel> &p,
const Point_2<Kernel> &q);

/*!
returns true iff `p` is lexicographically greater or equal to `q`.
\relates Point_2
*/
bool operator>=(const Point_2<Kernel> &p,
const Point_2<Kernel> &q);

/*!
returns the difference vector between `q` and `p`.
You can substitute `ORIGIN` for either `p` or `q`,
but not for both.
\relates Point_2
*/
Vector_2<Kernel> operator-(const Point_2<Kernel> &p,
const Point_2<Kernel> &q);


/// \ingroup Kernel_operator_plus
/// @{

/*!
returns the point obtained by translating `p` by the
vector `v`.
\relates Point_2
*/
Point_2<Kernel> operator+(const Point_2<Kernel> &p,
const Vector_2<Kernel> &v);

/// @}

/// \ingroup Kernel_operator_minus
/// @{

/*!
returns the point obtained by translating `p` by the
vector -`v`.
\relates Point_2
*/
Point_2<Kernel> operator-(const Point_2<Kernel> &p,
const Vector_2<Kernel> &v);

/// @}

} /* end namespace CGAL */
