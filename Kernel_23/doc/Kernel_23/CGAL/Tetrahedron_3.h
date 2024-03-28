namespace CGAL {

/*!
\ingroup kernel_classes3

An object `t` of the class `Tetrahedron_3` is an oriented
tetrahedron in the three-dimensional Euclidean space \f$ \E^3\f$.

It is defined by four vertices \f$ p_0\f$, \f$ p_1\f$, \f$ p_2\f$ and \f$ p_3\f$.
The orientation of a tetrahedron is the orientation of its four
vertices. That means it is positive when \f$ p_3\f$ is on the positive
side of the plane defined by \f$ p_0\f$, \f$ p_1\f$ and \f$ p_2\f$.

The tetrahedron itself splits the space \f$ \E^3\f$ in a <I>positive</I> and
a <I>negative</I> side.

The boundary of a tetrahedron splits the space in two open regions, a
bounded one and an unbounded one.

\cgalModels{Kernel::Tetrahedron_3}

*/
template< typename Kernel >
class Tetrahedron_3 {
public:

/// \name Creation
/// @{

/*!
introduces a tetrahedron `t` with vertices `p0`, `p1`, `p2` and `p3`.
\cgalEpicExact
*/
Tetrahedron_3(const Point_3<Kernel> &p0,
const Point_3<Kernel> &p1,
const Point_3<Kernel> &p2,
const Point_3<Kernel> &p3);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two tetrahedra `t` and `t2` are equal,
iff `t` and `t2` have the same orientation and
their sets (not sequences) of vertices are equal.
*/
bool operator==(const Tetrahedron_3<Kernel> &t2) const;

/*!
Test for inequality.
*/
bool operator!=(const Tetrahedron_3<Kernel> &t2) const;

/*!
returns the i'th vertex modulo 4 of `t`.
\cgalEpicExact
*/
Point_3<Kernel> vertex(int i) const;

/*!
returns `vertex(int i)`.
\cgalEpicExact
*/
Point_3<Kernel> operator[](int i) const;

/// @}

/// \name Predicates
/// @{

/*!
Tetrahedron `t` is degenerate, if the vertices are coplanar.
*/
bool is_degenerate() const;

/*!

*/
Orientation orientation() const;

/*!
\pre `t` is not degenerate.
*/
Oriented_side oriented_side(const Point_3<Kernel> &p) const;

/*!
\pre `t` is not degenerate.
*/
Bounded_side bounded_side(const Point_3<Kernel> &p) const;

/// @}

/// \name Convencience Boolean Functions
/// @{

/*!

*/
bool has_on_positive_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_negative_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_boundary(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_bounded_side(const Point_3<Kernel> &p) const;

/*!

*/
bool has_on_unbounded_side(const Point_3<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns the signed volume of `t`.
*/
Kernel::FT volume() const;

/*!
returns a bounding box containing `t`.
\cgalEpicExact
*/
Bbox_3 bbox() const;

/*!
returns the tetrahedron obtained by applying `at` on the three
vertices of `t`.
*/
Tetrahedron_3<Kernel> transform(const Aff_transformation_3<Kernel> &at) const;

/// @}

}; /* end Tetrahedron_3 */
} /* end namespace CGAL */
