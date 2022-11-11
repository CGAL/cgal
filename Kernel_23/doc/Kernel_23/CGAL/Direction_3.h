namespace CGAL {

/*!
\ingroup kernel_classes3

An object of the class `Direction_3` is a vector in the three-dimensional
vector space \f$ \mathbb{R}^3\f$ where we forget about their length. They can be
viewed as unit vectors, although there is no normalization internally,
since this is error prone. Directions are used whenever the length of
a vector does not matter.
They also characterize a set of parallel lines that have the same orientation
or the direction normal to parallel planes that have the same orientation.
For example, you can ask for the direction
orthogonal to an oriented plane, or the direction of an oriented line.

\cgalModels `Kernel::Direction_3`

*/
template< typename Kernel >
class Direction_3 {
public:

/// \name Creation
/// @{

/*!
introduces a direction `d` initialized with the
direction of vector `v`.
*/
Direction_3(const Vector_3<Kernel> &v);

/*!
introduces the direction `d` of line `l`.
*/
Direction_3(const Line_3<Kernel> &l);

/*!
introduces the direction `d` of ray `r`.
*/
Direction_3(const Ray_3<Kernel> &r);

/*!
introduces the direction `d` of segment `s`.
*/
Direction_3(const Segment_3<Kernel> &s);

/*!
introduces a direction `d` initialized with the direction
from the origin to the point with %Cartesian coordinates \f$ (x, y, z)\f$.
*/
Direction_3(const Kernel::RT &x, const Kernel::RT &y, const Kernel::RT &z);

/// @}

/// \name Operations
/// @{

/*!
returns values, such that `d``== Direction_3<Kernel>(delta(0),delta(1),delta(2))`.
\pre \f$ 0 \leq i \leq2\f$.
*/
Kernel::RT delta(int i) const;

/*!
returns `delta(0)`.
*/
Kernel::RT dx() const;

/*!
returns `delta(1)`.
*/
Kernel::RT dy() const;

/*!
returns `delta(2)`.
*/
Kernel::RT dz() const;

/*!
Test for equality.
*/
bool operator==(const Direction_3<Kernel> &e) const;

/*!
Test for inequality.
*/
bool operator!=(const Direction_3<Kernel> &e) const;

/*!
The direction opposite to `d`.
*/
Direction_3<Kernel> operator-() const;

/*!
returns a vector that has the same direction as `d`.
*/
Vector_3<Kernel> vector() const;

/*!
returns the direction obtained by applying `t` on `d`.
*/
Direction_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const;

/// @}

}; /* end Direction_3 */
} /* end namespace CGAL */
