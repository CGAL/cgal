namespace CGAL {
/*!
\ingroup kernel_classes2

An object `d` of the class `Direction_2` is a vector in the two-dimensional
vector space \f$ \mathbb{R}^2\f$ where we forget about its length. They can be
viewed as unit vectors, although there is no normalization internally,
since this is error prone. Directions are used whenever the length of
a vector does not matter.
They also characterize a set of parallel oriented lines that have the same
orientations.
For example, you can ask for the direction
orthogonal to an oriented plane, or the direction of an oriented line.
Further, they can be used to indicate angles. The slope of a direction
is `dy()`/`dx()`.

\cgalModels{Kernel::Direction_2}

*/
template< typename Kernel >
class Direction_2 {
public:

/// \name Creation
/// @{

/*!
introduces the direction `d` of vector `v`.
\cgalEpicExact
*/
Direction_2(const Vector_2<Kernel> &v);

/*!
introduces the direction `d` of line `l`.
*/
Direction_2(const Line_2<Kernel> &l);

/*!
introduces the direction `d` of ray `r`.
*/
Direction_2(const Ray_2<Kernel> &r);

/*!
introduces the direction `d` of segment `s`.
*/
Direction_2(const Segment_2<Kernel> &s);

/*!
introduces a direction `d` passing through the origin
and the point with %Cartesian coordinates \f$ (x, y)\f$.
\cgalEpicExact
*/
Direction_2(const Kernel::RT &x, const Kernel::RT &y);

/// @}

/// \name Operations
/// There is a total order on directions. We compare the angles
/// between the positive \f$ x\f$-axis and the directions in
/// counterclockwise order.
/// @{

/*!
returns values, such that `d``== Direction_2<Kernel>(delta(0),delta(1))`.
\pre `0 <= i <= 1`.

\cgalEpicExact
*/
Kernel::RT delta(int i) const;

/*!
returns `delta(0)`.
\cgalEpicExact
*/
Kernel::RT dx() const;

/*!
returns `delta(1)`.
\cgalEpicExact
*/
Kernel::RT dy() const;

/*!

*/
bool operator==(const Direction_2<Kernel> &e) const;

/*!

*/
bool operator!=(const Direction_2<Kernel> &e) const;

/*!

*/
bool operator<(const Direction_2<Kernel> &e) const;

/*!

*/
bool operator>(const Direction_2<Kernel> &e) const;

/*!

*/
bool operator<=(const Direction_2<Kernel> &e) const;

/*!

*/
bool operator>=(const Direction_2<Kernel> &e) const;

/*!
returns true, iff `d` is not equal to `d1`, and
while rotating counterclockwise starting at `d1`,
`d` is reached strictly before `d2` is reached.
Note that true is returned if `d1` == `d2`, unless
also `d` == `d1`.

*/
bool counterclockwise_in_between(const Direction_2<Kernel> &d1,
const Direction_2<Kernel> &d2) const;

/*!
The direction opposite to `d`.
*/
Direction_2<Kernel> operator-() const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns a vector that has the same direction as `d`.
*/
Vector_2<Kernel> vector() const;

/*!
returns the direction obtained by applying `t` on `d`.
*/
Direction_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const;

/// @}

}; /* end Direction_2 */
} /* end namespace CGAL */
