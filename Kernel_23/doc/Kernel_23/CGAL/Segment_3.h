namespace CGAL {

/*!
\ingroup kernel_classes3

An object `s` of the data type `Segment_3` is a directed
straight line segment in the three-dimensional Euclidean space \f$ \E^3\f$,
that is a straight line segment \f$ [p,q]\f$ connecting two points \f$ p,q \in
\R^3\f$. The segment is topologically closed, i.e.\ the end
points belong to it. Point `p` is called the <I>source</I> and `q`
is called the <I>target</I> of `s`. The length of `s` is the
Euclidean distance between `p` and `q`. Note that there is only a function
to compute the square of the length, because otherwise we had to
perform a square root operation which is not defined for all
number types, which is expensive, and may not be exact.

\cgalModels `Kernel::Segment_3`
\cgalModels `Hashable` if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`

*/
template< typename Kernel >
class Segment_3 {
public:

/// \name Creation
/// @{

/*!
introduces a segment `s` with source `p`
and target `q`. It is directed from the source towards
the target.
*/
Segment_3(const Point_3<Kernel> &p, const Point_3<Kernel> &q);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: Two segments are equal, iff their sources and
targets are equal.
*/
bool operator==(const Segment_3<Kernel> &q) const;

/*!
Test for inequality.
*/
bool operator!=(const Segment_3<Kernel> &q) const;

/*!
returns the source of `s`.
*/
Point_3<Kernel> source() const;

/*!
returns the target of `s`.
*/
Point_3<Kernel> target() const;

/*!
returns the point of `s` with smallest coordinate (lexicographically).
*/
Point_3<Kernel> min() const;

/*!
returns the point of `s` with largest coordinate (lexicographically).
*/
Point_3<Kernel> max() const;

/*!
returns source or target of `s`: `vertex(0)` returns
the source, `vertex(1)` returns the target.
The parameter `i` is taken modulo 2, which gives
easy access to the other vertex.
*/
Point_3<Kernel> vertex(int i) const;

/*!
returns `vertex(i)`.
*/
Point_3<Kernel> point(int i) const;

/*!
returns `vertex(i)`.
*/
Point_3<Kernel> operator[](int i) const;

/*!
returns the squared length of `s`.
*/
Kernel::FT squared_length() const;

/*!
returns the vector `s.target()` - `s`.`source()`.
*/
Vector_3<Kernel> to_vector() const;

/*!
returns the direction from source to target.
*/
Direction_3<Kernel> direction() const;

/*!
returns a segment with source and target interchanged.
*/
Segment_3<Kernel> opposite() const;

/*!
returns the line `l` passing through `s`. Line `l` has the
same orientation as segment `s`, that is
from the source to the target of `s`.
*/
Line_3<Kernel> supporting_line() const;

/*!
segment `s` is degenerate, if source and target fall together.
*/
bool is_degenerate() const;

/*!
A point is on `s`, iff it is equal to the source or target
of `s`, or if it is in the interior of `s`.
*/
bool has_on(const Point_3<Kernel> &p) const;

/*!
returns a bounding box containing `s`.
*/
Bbox_3 bbox() const;

/*!
returns the segment obtained by applying `t` on the source
and the target of `s`.
*/
Segment_3<Kernel> transform(const Aff_transformation_3<Kernel> &t) const;

/// @}

}; /* end Segment_3 */
} /* end namespace CGAL */
