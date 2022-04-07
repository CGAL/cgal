namespace CGAL {

/*!
\ingroup kernel_classes2

An object `r` of the data type `Ray_2` is a directed
straight ray in the two-dimensional Euclidean plane \f$ \E^2\f$. It starts
in a point called the <I>source</I> of `r` and goes to infinity.

\cgalModels `Kernel::Ray_2`

*/
template< typename Kernel >
class Ray_2 {
public:

/// \name Creation
/// @{

/*!
introduces a ray `r`
with source `p` and passing through point `q`.
*/
Ray_2(const Point_2<Kernel> &p, const Point_2<Kernel>&q);

/*!
introduces a ray `r` starting at source `p` with
direction `d`.
*/
Ray_2(const Point_2<Kernel> &p, const Direction_2<Kernel> &d);

/*!
introduces a ray `r` starting at source `p` with
the direction of `v`.
*/
Ray_2(const Point_2<Kernel> &p, const Vector_2<Kernel> &v);

/*!
introduces a ray `r` starting at source `p` with
the same direction as `l`.
*/
Ray_2(const Point_2<Kernel> &p, const Line_2<Kernel> &l);

/// @}

/// \name Operations
/// @{

/*!
Test for equality: two rays are equal, iff they have the same
source and the same direction.
*/
bool operator==(const Ray_2<Kernel> &h) const;

/*!
Test for inequality.
*/
bool operator!=(const Ray_2<Kernel> &h) const;

/*!
returns the source of `r`.
*/
Point_2<Kernel> source() const;

/*!
returns a point on `r`. `point(0)` is the source,
`point(i)`, with `i>0`, is different from the
source. \pre \f$ i \geq0\f$.
*/
Point_2<Kernel> point(const Kernel::FT i) const;

/*!
returns the direction of `r`.
*/
Direction_2<Kernel> direction() const;

/*!
returns a vector giving the direction of `r`.
*/
Vector_2<Kernel> to_vector() const;

/*!
returns the line supporting `r` which has the same direction.
*/
Line_2<Kernel> supporting_line() const;

/*!
returns the ray with the same source and the opposite direction.
*/
Ray_2<Kernel> opposite() const;

/// @}

/// \name Predicates
/// @{

/*!
ray `r` is degenerate, if the source and the second defining
point fall together (that is if the direction is degenerate).
*/
bool is_degenerate() const;

/*!

*/
bool is_horizontal() const;

/*!

*/
bool is_vertical() const;

/*!
A point is on `r`, iff it is equal to the source
of `r`, or if it is in the interior of `r`.
*/
bool has_on(const Point_2<Kernel> &p) const;

/*!
checks if point `p` is on `r`. This function is faster
than function `has_on()` if the precondition
checking is disabled.
\pre `p` is on the supporting line of `r`.
*/
bool collinear_has_on(const Point_2<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
returns the ray obtained by applying `t` on the source
and on the direction of `r`.
*/
Ray_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const;

/// @}

}; /* end Ray_2 */
} /* end namespace CGAL */
