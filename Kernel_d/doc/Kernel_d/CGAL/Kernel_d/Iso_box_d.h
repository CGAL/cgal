namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An object \f$ b\f$ of the data type `Iso_box_d` is an
iso-box in the Euclidean space \f$ \E^d\f$ with edges parallel to the
axes of the coordinate system.

*/
template< typename Kernel >
class Iso_box_d {
public:

/// \name Creation
/// @{

/*!
introduces an iso-oriented iso-box `b` with diagonal
opposite vertices \f$ p\f$ and \f$ q\f$.

*/
Iso_box_d(const Point_d<Kernel>& p,
const Point_d<Kernel> &q);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two iso-oriented cuboid are equal, iff their
lower left and their upper right vertices are equal.
*/
bool operator==(const Iso_box_d<Kernel>& b2) const;

/*!
Test for inequality.
*/
bool operator!=(const Iso_box_d<Kernel>& b2) const;

/*!
returns the smallest vertex of `b`.
*/
const Point_d<Kernel>& min() const;

/*!
returns the largest vertex of `b`.
*/
const Point_d<Kernel>& max() const;

/// @}

/// \name Predicates
/// @{

/*!
`b` is degenerate, if all vertices
are collinear.
*/
bool is_degenerate() const;

/*!
returns either `ON_UNBOUNDED_SIDE`,
`ON_BOUNDED_SIDE`, or the constant
`ON_BOUNDARY`,
depending on where point \f$ p\f$ is.
*/
Bounded_side bounded_side(const Point_d<Kernel>& p) const;

/*!

*/
bool has_on_boundary(const Point_d<Kernel>& p) const;

/*!

*/
bool has_on_bounded_side(const Point_d<Kernel>& p) const;

/*!

*/
bool has_on_unbounded_side(const Point_d<Kernel>& p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns the volume of `b`.
*/
Kernel_d::FT volume() const;

/// @}

}; /* end Iso_box_d */
} /* end namespace CGAL */
