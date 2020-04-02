namespace CGAL {
/*!
\ingroup kernel_classes3

An object `c` of the data type `Iso_cuboid_3` is a
cuboid in the Euclidean space \f$ \E^3\f$ with edges parallel to the \f$ x\f$,
\f$ y\f$ and \f$ z\f$ axis of the coordinate system.

Although they are represented in a canonical form by only two
vertices, namely the lexicographically smallest and largest vertex
with respect to %Cartesian \f$ xyz\f$ coordinates, we provide
functions for "accessing" the other vertices as well.

Iso-oriented cuboids and bounding boxes are quite similar. The
difference however is that bounding boxes have always double coordinates,
whereas the coordinate type of an iso-oriented cuboid is chosen by
the user.

\cgalModels `Kernel::IsoCuboid_3`
\cgalModels `Hashable` if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`

*/
template< typename Kernel >
class Iso_cuboid_3 {
public:

/// \name Creation
/// @{

/*!
introduces an iso-oriented cuboid `c` with diagonal
opposite vertices `p` and `q`. Note that the object is
brought in the canonical form.
*/
Iso_cuboid_3(const Point_3<Kernel> &p,
const Point_3<Kernel> &q);

/*!
introduces an iso-oriented cuboid `c` with diagonal
opposite vertices `p` and `q`. The `int` argument value
is only used to distinguish the two overloaded functions.
\pre `p.x()<=q.x()`, `p.y()<=q.y()`and `p.z()<=q.z()`.
*/
Iso_cuboid_3(const Point_3<Kernel> &p,
const Point_3<Kernel> &q, int);

/*!
introduces an iso-oriented cuboid `c` whose
minimal \f$ x\f$ coordinate is the one of `left`, the
maximal \f$ x\f$ coordinate is the one of `right`, the
minimal \f$ y\f$ coordinate is the one of `bottom`, the
maximal \f$ y\f$ coordinate is the one of `top`, the
minimal \f$ z\f$ coordinate is the one of `far`, the
maximal \f$ z\f$ coordinate is the one of `close`.
*/
Iso_cuboid_3(const Point_3<Kernel> &left,
const Point_3<Kernel> &right,
const Point_3<Kernel> &bottom,
const Point_3<Kernel> &top,
const Point_3<Kernel> &far,
const Point_3<Kernel> &close);

/*!
introduces an iso-oriented cuboid `c` with diagonal
opposite vertices
(`min_hx/hw`, `min_hy/hw`, `min_hz/hw`) and
(`max_hx/hw`, `max_hy/hw`, `max_hz/hw`).
\pre `hw` \f$ \neq\f$ 0.
*/
Iso_cuboid_3(
const Kernel::RT& min_hx, const Kernel::RT& min_hy, const Kernel::RT& min_hz,
const Kernel::RT& max_hx, const Kernel::RT& max_hy, const Kernel::RT& max_hz,
const Kernel::RT& hw = RT(1));

/*!
If `Kernel::RT` is constructible from double,
introduces an iso-oriented cuboid from `bbox`.
*/
Iso_cuboid_3(const Bbox_3& bbox);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two iso-oriented cuboid are equal, iff their
lower left and their upper right vertices are equal.
*/
bool operator==(const Iso_cuboid_3<Kernel> &c2) const;

/*!
Test for inequality.
*/
bool operator!=(const Iso_cuboid_3<Kernel> &c2) const;

/*!
returns the i'th vertex modulo 8 of `c`.
starting with the lower left vertex.
*/
Point_3<Kernel> vertex(int i) const;

/*!
returns `vertex(i)`, as indicated in the figure below:
\image html IsoCuboid.png
\image latex IsoCuboid.png
*/
Point_3<Kernel> operator[](int i) const;

/*!
returns the smallest vertex of `c` (= `vertex(0)`).
*/
Point_3<Kernel> min() const;

/*!
returns the largest vertex of `c` (= `vertex(7)`).
*/
Point_3<Kernel> max() const;

/*!
returns smallest %Cartesian
\f$ x\f$-coordinate in `c`.
*/
Kernel::FT xmin() const;

/*!
returns smallest %Cartesian
\f$ y\f$-coordinate in `c`.
*/
Kernel::FT ymin() const;

/*!
returns smallest %Cartesian
\f$ z\f$-coordinate in `c`.
*/
Kernel::FT zmin() const;

/*!
returns largest %Cartesian
\f$ x\f$-coordinate in `c`.
*/
Kernel::FT xmax() const;

/*!
returns largest %Cartesian
\f$ y\f$-coordinate in `c`.
*/
Kernel::FT ymax() const;

/*!
returns largest %Cartesian
\f$ z\f$-coordinate in `c`.
*/
Kernel::FT zmax() const;

/*!
returns `i`-th %Cartesian coordinate of
the smallest vertex of `c`.
\pre \f$ 0 \leq i \leq2\f$.
*/
Kernel::FT min_coord(int i) const;

/*!
returns `i`-th %Cartesian coordinate of
the largest vertex of `c`.
\pre \f$ 0 \leq i \leq2\f$.
*/
Kernel::FT max_coord(int i) const;

/// @}

/// \name Predicates
/// @{

/*!
`c` is degenerate, if all vertices
are coplanar.
*/
bool is_degenerate() const;

/*!
returns either \ref ON_UNBOUNDED_SIDE,
\ref ON_BOUNDED_SIDE, or the constant
\ref ON_BOUNDARY,
depending on where point `p` is.
*/
Bounded_side bounded_side(const Point_3<Kernel> &p) const;

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
returns the volume of `c`.
*/
Kernel::FT volume() const;

/*!
returns a bounding box containing `c`.
*/
Bbox_3 bbox() const;

/*!
returns the iso-oriented cuboid obtained by applying `t` on
the smallest and the largest of `c`.
\pre The angle at a rotation must be a multiple of \f$ \pi/2\f$, otherwise the resulting cuboid does not have the same size. Note that rotating about an arbitrary angle can even result in a degenerate iso-oriented cuboid.
*/
Iso_cuboid_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const;

/// @}

}; /* end Iso_cuboid_3 */
} /* end namespace CGAL */
