namespace CGAL {

/*!
\ingroup kernel_classes2

An object `t` of the class `Triangle_2` is a triangle
in the two-dimensional Euclidean plane \f$ \E^2\f$.
Triangle `t` is oriented, i.e., its boundary has
clockwise or counterclockwise orientation. We call the side to the left
of the boundary the positive side and the side to the right of the
boundary the negative side.

The boundary of a triangle splits the plane in
two open regions, a bounded one and an unbounded one.

\cgalModels `Kernel::Triangle_2`

*/
template< typename Kernel >
class Triangle_2 {
public:

/// \name Creation
/// @{

/*!
introduces a triangle `t` with vertices `p`, `q` and `r`.
*/
Triangle_2(const Point_2<Kernel> &p,
const Point_2<Kernel> &q,
const Point_2<Kernel> &r);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two triangles are equal, iff there exists a
cyclic permutation of the vertices of \f$ t2\f$, such that they are
equal to the vertices of `t`.
*/
bool operator==(const Triangle_2<Kernel> &t2) const;

/*!
Test for inequality.
*/
bool operator!=(const Triangle_2<Kernel> &t2) const;

/*!
returns the i'th vertex modulo 3 of `t`.
*/
Point_2<Kernel> vertex(int i) const;

/*!
returns `vertex(i)`.
*/
Point_2<Kernel> operator[](int i) const;

/// @}

/// \name Predicates
/// For convenience we provide the following Boolean functions:
/// @{

/*!
triangle `t` is degenerate, if the vertices are collinear.
*/
bool is_degenerate() const;

/*!
returns the orientation of `t`.
*/
Orientation orientation() const;

/*!
returns
`ON_ORIENTED_BOUNDARY`, or
`POSITIVE_SIDE`,
or the constant
`ON_NEGATIVE_SIDE`,
determined by the position of point `p`.
\pre `t` is not degenerate.
*/
Oriented_side oriented_side(const Point_2<Kernel> &p) const;

/*!
returns the constant `ON_BOUNDARY`,
`ON_BOUNDED_SIDE`, or else
`ON_UNBOUNDED_SIDE`,
depending on where point `p` is.
\pre `t` is not degenerate.
*/
Bounded_side bounded_side(const Point_2<Kernel> &p) const;

/*!

*/
bool has_on_positive_side(const Point_2<Kernel> &p) const;

/*!

*/
bool has_on_negative_side(const Point_2<Kernel> &p) const;

/*!

*/
bool has_on_boundary(const Point_2<Kernel> &p) const;

/*!

*/
bool has_on_bounded_side(const Point_2<Kernel> &p) const;

/*!
\pre `t` is not degenerate.
*/
bool has_on_unbounded_side(const Point_2<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns a triangle where the boundary is oriented the other
way round (this flips the positive and the negative side, but
not the bounded and unbounded side).
*/
Triangle_2<Kernel> opposite();

/*!
returns the signed area of `t`.
*/
Kernel::FT area() const;

/*!
returns a bounding box containing `t`.
*/
Bbox_2 bbox() const;

/*!
returns the triangle obtained by applying \f$ at\f$ on the three
vertices of `t`.
*/
Triangle_2<Kernel> transform(const Aff_transformation_2<Kernel> &at) const;

/// @}

}; /* end Triangle_2 */
} /* end namespace CGAL */
