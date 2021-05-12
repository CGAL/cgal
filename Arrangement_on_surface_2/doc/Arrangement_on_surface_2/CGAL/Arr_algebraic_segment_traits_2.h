
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

The traits class `Arr_algebraic_segment_traits_2` is a model of the `ArrangementTraits_2`
concept that handles planar algebraic curves of arbitrary degree,
and \f$ x\f$-monotone of such curves.
A planar (real) <I>algebraic curve</I>
is the vanishing set of a polynomial in two variables, that is, the
curve is defined by the defining equation
\f[ f(x):=\sum_{i+j\leq n} a_{ij} x^i y^j =0, \f]
where \f$ n\f$ is the degree of the curve.

The traits class allows the construction of algebraic curves,
by specifying their implicit equation. \f$ x\f$-monotone and vertical segments
of a curve can also be defined; unbounded curves and segments are supported.
The template parameter `Coefficient` defines
the innermost coefficient type of the polynomials. Currently,
the types `leda::integer` and `CORE::BigInt` are supported as well
as any instance of `CGAL::Sqrt_extension` that is instantiated with
one of the integral types above.

\cgalModels `ArrangementTraits_2`


*/
template< typename Coefficient >
class Arr_algebraic_segment_traits_2 {
public:

/// \name Types
/// @{

/*!
Value to specify whether a point should be in the interior
of a segment, or its minimal point,
or its maximal point in lexicographic order.
*/
enum Site_of_point { POINT_IN_INTERIOR = 0, MIN_ENDPOINT = -1, MAX_ENDPOINT = 1 };

/*!
the type for bivariate polynomials,
with innermost coefficient type `Coefficient`.
Constitutes a model of the concept `Polynomial_d`
with two variables.

\sa `CGAL::Polynomial_d`
*/
typedef unspecified_type Polynomial_2;

/*!
model for the concept
`AlgebraicKernel_1`
*/
typedef unspecified_type Algebraic_kernel_1;

/*!
represents coordinates of points.
Typedef from `Algebraic_kernel_1::Algebraic_real_1`
*/
typedef unspecified_type Algebraic_real_1;

/*!
Typedef from `Algebraic_kernel_1::Bound`
*/
typedef unspecified_type Bound;

/// @}

/// \name Accessing functor objects
/// @{

/*!

*/
Construct_curve_2 construct_curve_2_object() const;

/*!

*/
Construct_point_2 construct_point_2_object() const;

/*!

*/
Construct_x_monotone_segment_2 construct_x_monotone_segment_2_object() const;

/// @}

/*!


*/
class Construct_curve_2 {
public:

/// \name Object Creation Functors
/// @{

/*!
Returns a `Curve_2` object that represents the curve defined by
the polynomial `p`
*/
Curve_2 operator() (Polynomial_2 p);

/*!
Returns a `Curve_2` object specified by `s`.
The passed string represents the defining polynomial of the curve
and must be given in a MAPLE-readable format using "x" as first
and "y" as second variable, e.g.,
"(x^3*y-2*x)*(-6*x-y^3*x^6)"
for integer coefficients, and "3/2*x*y^4-5/7*x^2+3/1"
for rational coefficients.
*/
Curve_2 operator() (std::string s);

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_curve_2 */

/*!


*/
class Construct_point_2 {
public:

/// \name Object Creation Functors
/// @{

/*!
Returns a `Point_2` object that represents the `arcno`-th
point in the fiber of `cv` at \f$ x\f$-coordinate `x`,
counted from the bottom, starting with zero.
\pre (`cv` must not have a vertical line at `x`,
and \f$ 0\leq arcno < c\f$, where \f$ c\f$ is the number of points
in the fiber of `cv` at `x`.)
*/
Point_2 operator() (Algebraic_real_1 x, Curve_2 cv, int arcno);

/*!
Returns a `Point_2` object that represents the
point on `xcv` at \f$ x\f$-coordinate `x`
\pre (`x` is in the \f$ x\f$-range of `xcv`.)
*/
Point_2 operator() (Algebraic_real_1 x, X_monotone_curve_2 xcv);

/*!
Returns a `Point_2` object that represents (x,y)
*/
Point_2 operator() (Algebraic_real_1 x, Algebraic_real_1 y);

/*!
Returns a `Point_2` object that represents (x,y)
*/
Point_2 operator() (Coefficient x, Coefficient y);

/*!
Returns a `Point_2` object that represents (x,y)
*/
Point_2 operator() (Bound x, Bound y);

/*!
Returns a `Point_2` object that represents (x,y)
*/
Point_2 operator() (int x, int y);

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_point_2 */

/*!


*/
class Construct_x_monotone_segment_2 {
public:

/// \name Object Creation Functors
/// @{

/*!
Writes a sequence of `X_monotone_curve_2` objects (terminal
segments) into `out`. These terminal segments compose an
\f$ x\f$-monotone (or vertical) segment of the curve `cv` that
starts in `end_min`, and ends in `end_max`.
\pre (`end_min` must have a unique \f$ x\f$-monotone segment
to its right, or `end_max` must have a unique
\f$ x\f$-monotone segment to its left. Furthermore,
`end_min` and `end_max` must be connected
by an \f$ x\f$-monotone segment of `cv`)
*/
template<class OutputIterator> OutputIterator
operator() (Curve_2 cv, Point_2 end_min, Point_2 end_max,
OutputIterator out);

/*!
Writes a sequence of `X_monotone_curve_2` objects into `out`.
These segments form an \f$ x\f$-monotone (or vertical)
segment of the curve `cv`.

If `site_of_p==POINT_IN_INTERIOR`, the maximal segment is
returned that contains `p` in its interior.

returned that contains `p` as its left endpoint.

returned that contains `p` as its left endpoint.
\pre (If `site_of_p==POINT_IN_INTERIOR`, `p`
must be an interior point of an \f$ x\f$-monotone or a vertical
segment.
must either have a unique \f$ x\f$-monotone segment to the right,
or a vertical segment from `p` upwards.
must either have a unique \f$ x\f$-monotone segment to the left,
or a vertical segment from `p` downwards.)
*/
template<class OutputIterator> OutputIterator
operator() (Curve_2 cv, Point_2 p,
Site_of_point site_of_p,
OutputIterator out);

/*!
Writes a sequence of `X_monotone_curve_2` objects into `out`.
These segments form a straight-line segment connecting
the points `p` and `q`. If `p` and `q` share the
same \f$ x\f$-coordinate, the constructed vertical segment consists of
only one `X_monotone_curve_2` object and can be computed
efficiently. In the non-vertical case,
the construction is only possible if `p` and `q`
have both rational x- and y-coordinates.
\pre (`p` must not be equal to `q`.)

*/
template<class OutputIterator> OutputIterator
operator() (Point_2 p, Point_2 q,
OutputIterator out);

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_x_monotone_segment_2 */

/*!


Models the `ArrangementTraits_2::Curve_2` concept.
Represents algebraic curves. Internally, the type stores
topological-geometric information about the particular curve.
In order to use internal caching, instances should only be created
using the `Construct_curve_2` functor of the traits class.

*/
class Curve_2 {
public:

/// \name Modifiers
/// @{

/*!
returns the defining polynomial of the curve.
*/
Polynomial_2 polynomial () const;

/// @}

}; /* end Arr_algebraic_segment_traits_2::Curve_2 */

/*!


Models the `ArrangementBasicTraits_2::Point_2` concept.
Represents points in \f$ \mathbb{R}^2\f$. Intersection points of algebraic curves
are in general non-rational, so we need a data structure that is capable
of representing arbitrary points with algebraic coordinates.

The traits class represents algebraic coordinates by the type
`Algebraic_real_1`, which is a model of the `AlgebraicReal_1`
concept.
A point \f$ p\f$ is stored by a triple \f$ (x,cv,arcno)\f$,
where \f$ x\f$ is the \f$ x\f$-coordinate of a point, \f$ cv\f$ is an instance
of `Curve_2` that contains the point, (and has no vertical line at \f$ x\f$),
and \f$ arcno\f$ is an `int`, denoting that \f$ p\f$ is met as the \f$ arcno\f$-th
point when shooting a vertical ray at \f$ x\f$, starting from \f$ -\infty\f$
(where counting starts with \f$ 0\f$).

In addition to the methods listed below, the copy constructor and assignment
operator for `Point_2` objects are also supported.

The functor `Construct_point_2` constructs `Point_2` instances.

*/
class Point_2 {
public:

/// \name Modifiers
/// @{

/*!
returns the \f$ x\f$-coordinate of `p`.
*/
Algebraic_real_1 x () const;

/*!
returns the \f$ y\f$-coordinates of `p`.

<B>Attention:</B> As described above, points are not stored
by their \f$ y\f$-coordinate in `Algebraic_real_1` representation. In fact,
this representation must be computed on demand, and might become quite
costly for points defined by high-degree polynomials. Therefore, it is
recommended to avoid to call this function as much as possible.
*/
Algebraic_real_1 y () const;

/*!
returns a `Curve_2` instance that `p`is part of.
*/
Curve_2 curve () const;

/*!
returns the arc number of `p`.
*/
int arcno () const;

/*!
returns double-approximations of the \f$ x\f$- and \f$ y\f$-coordinates.
*/
std::pair<double,double> to_double () const;

/// @}

}; /* end Arr_algebraic_segment_traits_2::Point_2 */

/*!


Models the `ArrangementBasicTraits_2::X_monotone_curve_2` concept.
Represents terminal segments of an algebraic curves,
that means vertical segments or \f$ x\f$-monotone segments with no critical
\f$ x\f$-coordinate in the interior of their \f$ x\f$-range.
Terminal segments might either be bounded or unbounded.
By definition, each interior point of
a non-vertical segment has the same arc number (see the documentation of
type `Point_2` above, which is called the <I>arc number</I> of the segment
(note the arc number at the endpoints might differ).
Such segments are represented internally by a 4-tuple \f$ (p,q,cv,arcno)\f$,
where \f$ p\f$ and \f$ q\f$ are the endpoints, \f$ cv\f$ is the <I>supporting curve</I>
that the segment belongs to, and arcno is the arc number of the segment.

Arbitrary (weakly) \f$ x\f$-monotone segments are presented by a range
of `X_monotone_curve_2` instances, whose union equals the segment.
The functor `Construct_x_monotone_segment_2` allows their construction.
To construct all (maximal) terminal segments of a curve,
use the `Make_x_monotone_2` functor supplied by the traits class.

*/
class X_monotone_curve_2 {
public:

/// \name Modifiers
/// @{

/*!
returns the supporting algebraic curve of `s`.
*/
Curve_2 curve () const;

/*!
returns whether `s` is a vertical segment.
*/
bool is_vertical () const;

/*!
returns whether `s` has a finite endpoint on the left
*/
bool is_finite (CGAL::Arr_curve_end ce) const;

/*!
\pre (The corresponding curve end is finite)
*/
Point_2 curve_end (CGAL::Arr_curve_end ce) const;

/*!
returns the arc number of the segment.
\pre (The segment is non-vertical)
*/
int arcno () const;

/*!
returns the \f$ x\f$-coordinate of a vertical segment.
\pre (The segment is vertical)
*/
Algebraic_real_1 x () const;

/// @}

}; /* end Arr_algebraic_segment_traits_2::X_monotone_curve_2 */




}; /* end Arr_algebraic_segment_traits_2 */
} /* end namespace CGAL */
