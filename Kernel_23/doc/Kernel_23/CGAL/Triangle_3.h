namespace CGAL {

/*!
\ingroup kernel_classes3

An object `t` of the class `Triangle_3` is a triangle in
the three-dimensional Euclidean space \f$ \E^3\f$. As the triangle is not
a full-dimensional object there is only a test whether a point lies on
the triangle or not.

\cgalModels{Kernel::Triangle_3}

*/
template< typename Kernel >
class Triangle_3 {
public:

/// \name Creation
/// @{

/*!
introduces a triangle `t` with vertices `p`, `q` and `r`.
\cgalEpicExact
*/
Triangle_3(const Point_3<Kernel> &p,
const Point_3<Kernel> &q,
const Point_3<Kernel> &r);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two triangles `t1` and `t2` are equal, iff there
exists a cyclic permutation of the vertices of `t2`, such that
they are equal to the vertices of `t1`.
*/
bool operator==(const Triangle_3<Kernel> &t2) const;

/*!
Test for inequality.
*/
bool operator!=(const Triangle_3<Kernel> &t2) const;

/*!
returns the i-th vertex modulo 3 of `t`.
\cgalEpicExact
*/
Point_3<Kernel> vertex(int i) const;

/*!
returns `vertex(int i)`.
\cgalEpicExact
*/
Point_3<Kernel> operator[](int i) const;

/*!
returns the supporting plane of `t`, with same
orientation.
*/
Plane_3<Kernel> supporting_plane();

/// @}

/// \name Predicates
/// @{

/*!
`t` is degenerate if its vertices are collinear.
*/
bool is_degenerate() const;

/*!
A point is on `t`, if it is on a vertex, an edge or the
face of `t`.
*/
bool has_on(const Point_3<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns a square of the area of `t`.
*/
Kernel::FT squared_area() const;

/*!
returns a bounding box containing `t`.
\cgalEpicExact
*/
Bbox_3 bbox() const;

/*!
returns the triangle obtained by applying `at` on the three
vertices of `t`.
*/
Triangle_3<Kernel> transform(const Aff_transformation_3<Kernel> &at) const;

/// @}

}; /* end Triangle_3 */
} /* end namespace CGAL */
