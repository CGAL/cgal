namespace CGAL {

/*!
\ingroup kernel_classes3

An object of the class `Point_3` is a point in the three-dimensional 
Euclidean space \f$ \E^3\f$. 

Remember that `Kernel::RT` and `Kernel::FT` denote a 
RingNumberType and a FieldNumberType, respectively. For the kernel 
model `Cartesian<T>`, the two types are the same. For the 
kernel model `Homogeneous<T>`, `Kernel::RT` is equal 
to `T`, and `Kernel::FT` is equal to 
`Quotient<T>`. 

\cgalHeading{Operators}

The following operations can be applied on points: 

\sa `Kernel::Point_3` 

*/
template< typename Kernel >
class Point_3 {
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
introduces a point with %Cartesian coordinates\f$ (0,0,0)\f$. 
*/ 
Point_3(const Origin &ORIGIN); 

/*!
introduces a point `p` initialized to `(x,y,z)`. 
*/ 
Point_3(int x, int y, int z); 

/*!
introduces a point `p` initialized to `(x,y,z)`
provided `RT` supports it. 
*/ 
Point_3(double x, double y, double z); 

/*!
introduces a point `p` initialized to `(hx/hw,hy/hw, hz/hw)`. 
\pre `hw` \f$ \neq\f$ 0. 
*/ 
Point_3(const Kernel::RT &hx, const Kernel::RT &hy, const Kernel::RT &hz, const Kernel::RT &hw = RT(1)); 

/*!
introduces a point `p` initialized to `(x,y,z)`. 
*/ 
Point_3(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z); 

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: Two points are equal, iff their \f$ x\f$, \f$ y\f$ and \f$ z\f$ 
coordinates are equal. 
*/ 
bool operator==(const Point_3<Kernel> &q) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Point_3<Kernel> &q) const; 

/// @}

/// \name Coordinate Access
/// There are two sets of coordinate access functions, namely to the
/// homogeneous and to the %Cartesian coordinates. They can be used
/// independently from the chosen kernel model. Note that you do not
/// loose information with the homogeneous representation, because the
/// FieldNumberType is a quotient.
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
returns the %Cartesian \f$ x\f$ coordinate, that is `hx()`/`hw()`. 
*/ 
Kernel::FT x() const; 

/*!
returns the %Cartesian \f$ y\f$ coordinate, that is `hy()`/`hw()`. 
*/ 
Kernel::FT y() const; 

/*!
returns the %Cartesian \f$ z\f$ coordinate, that is `hz()`/`hw()`. 
*/ 
Kernel::FT z() const; 

/// @}

/// \name Convenience Operators
/// The following operations are for convenience and for compatibility
/// with code for higher dimensional points. Again they come in a
/// %Cartesian and in a homogeneous flavor.
/// @{

/*!
returns the i'th homogeneous coordinate of `p`, starting with 0. 
\pre \f$ 0\leq i \leq3\f$. 
*/ 
Kernel::RT homogeneous(int i) const; 

/*!
returns the i'th %Cartesian coordinate of `p`, starting with 0. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
Kernel::FT cartesian(int i) const; 

/*!
returns `cartesian(i)`. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
Kernel::FT operator[](int i) const; 

/*!
returns an iterator to the %Cartesian coordinates 
of `p`, starting with the 0th coordinate. 
*/ 
Cartesian_const_iterator cartesian_begin() const; 

/*!
returns an off the end iterator to the %Cartesian 
coordinates of `p`. 
*/ 
Cartesian_const_iterator cartesian_end() const; 

/*!
returns the dimension (the constant 3). 
*/ 
int dimension() const; 

/*!
returns a bounding box containing `p`. 
*/ 
Bbox_3 bbox() const; 

/*!
returns the point obtained by applying `t` on `p`. 
*/ 
Point_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const; 

/// @}

}; /* end Point_3 */

/*!
returns true iff `p` is lexicographically smaller than `q` 
(the lexicographical order being defined on the %Cartesian 
coordinates). 
\relates Point_3 
*/ 
bool operator<(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q); 

/*!
returns true iff `p` is lexicographically greater than `q`. 
\relates Point_3 
*/ 
bool operator>(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q); 

/*!
returns true iff `p` is lexicographically smaller or equal to 
`q`. 
\relates Point_3 
*/ 
bool operator<=(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q); 

/*!
returns true iff `p` is lexicographically greater or equal to 
`q`. 
\relates Point_3 
*/ 
bool operator>=(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q); 

/*!
returns the difference vector between `q` and `p`. 
You can substitute `ORIGIN` for either `p` or `q`, 
but not for both. 
\relates Point_3 
*/ 
Vector_3<Kernel> operator-(const Point_3<Kernel> &p, 
const Point_3<Kernel> &q); 

/// \ingroup Kernel_operator_plus

///@{

/*!
returns the point obtained by translating `p` by the 
vector `v`. 
\relates Point_3 
*/ 
Point_3<Kernel> operator+(const Point_3<Kernel> &p, 
const Vector_3<Kernel> &v); 

/// @}

/// \ingroup Kernel_operator_minus

///@{

/*!
returns the point obtained by translating `p` by the 
vector -`v`. 
\relates Point_3 
*/ 
Point_3<Kernel> operator-(const Point_3<Kernel> &p, 
const Vector_3<Kernel> &v); 

/// @}

} /* end namespace CGAL */
