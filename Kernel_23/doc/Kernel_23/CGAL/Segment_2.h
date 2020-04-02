namespace CGAL {

/*!
\ingroup kernel_classes2

An object `s` of the data type `Segment_2` is a directed
straight line segment in the two-dimensional Euclidean plane \f$ \E^2\f$, i.e.\ a
straight line segment \f$ [p,q]\f$ connecting two points \f$ p,q \in \mathbb{R}^2\f$.
The segment is topologically closed, i.e.\ the end
points belong to it. Point `p` is called the <I>source</I> and `q`
is called the <I>target</I> of `s`. The length of `s` is the
Euclidean distance between `p` and `q`. Note that there is only a function
to compute the square of the length, because otherwise we had to
perform a square root operation which is not defined for all
number types, which is expensive, and may not be exact.

\cgalModels `Kernel::Segment_2`
\cgalModels `Hashable` if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`

*/
template< typename Kernel >
class Segment_2 {
public:

/// \name Creation
/// @{

/*!
introduces a segment `s` with source `p`
and target `q`. The segment is directed from the source towards
the target.
*/
Segment_2(const Point_2<Kernel> &p, const Point_2<Kernel> &q);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: Two segments are equal, iff their sources and
targets are equal.
*/
bool operator==(const Segment_2<Kernel> &q) const;

/*!
Test for inequality.
*/
bool operator!=(const Segment_2<Kernel> &q) const;

/*!
returns the source of `s`.
*/
Point_2<Kernel> source() const;

/*!
returns the target of `s`.
*/
Point_2<Kernel> target() const;

/*!
returns the point of `s` with lexicographically smallest coordinate.
*/
Point_2<Kernel> min() const;

/*!
returns the point of `s` with lexicographically largest coordinate.
*/
Point_2<Kernel> max() const;

/*!
returns source or target of `s`: `vertex(0)` returns
the source of `s`, `vertex(1)` returns the target of `s`.
The parameter `i` is taken modulo 2, which gives
easy access to the other vertex.
*/
Point_2<Kernel> vertex(int i) const;

/*!
returns `vertex(i)`.
*/
Point_2<Kernel> point(int i) const;

/*!
returns `vertex(i)`.
*/
Point_2<Kernel> operator[](int i) const;

/*!
returns the squared length of `s`.
*/
Kernel::FT squared_length() const;

/*!
returns the direction from source to target of `s`.
*/
Direction_2<Kernel> direction() const;

/*!
returns the vector `s.target()` - `s.source()`.
*/
Vector_2<Kernel> to_vector() const;

/*!
returns a segment with source and target point interchanged.
*/
Segment_2<Kernel> opposite() const;

/*!
returns the line `l` passing through `s`. Line `l` has the
same orientation as segment `s`.
*/
Line_2<Kernel> supporting_line() const;

/// @}

/// \name Predicates
/// @{

/*!
segment `s` is degenerate, if source and target are equal.
*/
bool is_degenerate() const;

/*!

*/
bool is_horizontal() const;

/*!

*/
bool is_vertical() const;

/*!
A point is on `s`, iff it is equal to the source or target
of `s`, or if it is in the interior of `s`.
*/
bool has_on(const Point_2<Kernel> &p) const;

/*!
checks if point `p` is on segment `s`. This function is faster
than function `has_on()`.
\pre `p` is on the supporting line of `s`.
*/
bool collinear_has_on(const Point_2<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns a bounding box containing `s`.
*/
Bbox_2 bbox() const;

/*!
returns the segment obtained by applying `t` on the source
and the target of `s`.
*/
Segment_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const;

/// @}

}; /* end Segment_2 */
} /* end namespace CGAL */
