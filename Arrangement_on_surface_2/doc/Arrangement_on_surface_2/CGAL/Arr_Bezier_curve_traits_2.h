
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

The traits class `Arr_Bezier_curve_traits_2` is a model of the `ArrangementTraits_2`
concept that handles planar B&eacute;zier curves. A planar <I>B&eacute;zier curve</I>
\f$ B\f$ is a parametric curve defined by a sequence of <I>control points</I>
\f$ p_0, \ldots, p_n\f$ as follows:


\f{eqnarray*}{
B(t) = \left(X(t), Y(t)\right)
= \ccSum{k=0}{n}{p_k \cdot \frac{n!}{k! (n-k)!} \cdot
t^k (1-t)^{n-k}}\ .
\f}
where \f$ t \in [0, 1]\f$. The degree of the curve is therefore \f$ n\f$ -
namely, \f$ X(t)\f$ and \f$ Y(t)\f$ are polynomials of degree \f$ n\f$. B&eacute;zier curves
have numerous applications in computer graphics and solid modelling. They
are used, for example, in free-form sketches and for defining the true-type
fonts.

In our representation, we assume that the coordinates of all control
points are rational numbers (namely they are given as objects of the
`RatKernel::Point_2` type), so both \f$ X(t)\f$ and \f$ Y(t)\f$ are polynomials
with rational coefficients. The intersection points between curves are
however algebraic numbers, and their exact computation is time-consuming.
The traits class therefore contains a layer of geometric filtering that
performs all computation in an approximate manner whenever possible, and
it resorts to exact computations only when the approximate computation
fails to produce an unambiguous result.

We therefore require separate representations of the control points and
the intersection points. The `NtTraits` should be instantiated with a class
that defines nested `Integer`, `Rational` and `Algebraic` number
types and supports various operations on them, yielding certified computation
results (for example, in can convert rational numbers to algebraic numbers
and can compute roots of polynomials with integer coefficients).
The other template parameters, `RatKernel` and `AlgKernel` should be
geometric kernels templated with the `NtTraits::Rational` and
`NtTraits::Algebraic` number types, respectively. It is recommended to
instantiate the `CORE_algebraic_number_traits` class as the `NtTraits`
parameter, with `Cartesian<NtTraits::Rational>` and
`Cartesian<NtTraits::Algebraic>` instantiating the two kernel types,
respectively. The number types in this case are provided by the \core
library, with its ability to exactly represent simple algebraic numbers.

While `Arr_Bezier_curve_traits_2` models the concept
`ArrangementDirectionalXMonotoneTraits_2`, the implementation of
the `Are_mergeable_2` operation does not enforce the input curves
to have the same direction as a precondition. Moreover, `Arr_Bezier_curve_traits_2`
supports the merging of curves of opposite directions.

\cgalModels{ArrangementTraits_2,ArrangementDirectionalXMonotoneTraits_2}


*/
template< typename RatKernel, typename AlgKernel, typename NtTraits >
class Arr_Bezier_curve_traits_2 {
public:

/// \name Types
/// @{

/*!
the `NtTraits::Rational` type
(and also the `RatKernel::FT` type).
*/
typedef unspecified_type Rational;

/*!
the `NtTraits::Algebraic` type
(and also the `AlgKernel::FT` type).
*/
typedef unspecified_type Algebraic;

/// @}


/*!


The `Curve_2` class nested within the B&eacute;zier traits class is used
to represent a B&eacute;zier curve of arbitrary degree, which is defined by a
sequence of rational control points. In addition to the methods listed
below, the I/O operators \link PkgArrangementOnSurface2op_left_shift `operator<<` \endlink and \link PkgArrangementOnSurface2op_right_shift `operator>>` \endlink for
standard output-streams are also supported. The copy constructor and
assignment operator are supported as well.

*/
class Curve_2 {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Curve_2 ();

/*!
constructs a B&eacute;zier curve as defined by the given range of control
points. The value-type of `InputIterator` is `RatKernel::Point_2`.
\pre The input range must contain at least two control points.

*/
template <class InputIterator>
Curve_2 (InputIterator pts_begin, InputIterator pts_end);

/// @}

/// \name Access Functions
/// @{

/*!
returns the number of control points that define `B`.
*/
size_t number_of_control_points () const;

/*!
returns the \f$ k\f$th control point. Note that the first control point equals
the curve source, while the last control point equals its target. The rest
of the control points do not lie on the curve.
\pre \f$ k\f$ is smaller than the number of control points.
*/
typename RatKernel::Point_2 control_point (size_t k) const;

/*!
returns the point \f$ B(t)\f$ on the curve that corresponds to the given
rational parameter value.
*/
typename RatKernel::Point_2 operator() (const Rational& t) const;

/*!
returns the point \f$ B(t)\f$ on the curve that corresponds to the given
algebraic parameter value.
*/
typename AlgKernel::Point_2 operator() (const Algebraic& t) const;

/// @}

}; /* end Arr_Bezier_curve_traits_2::Curve_2 */


/*!

The `Point_2` class nested within the B&eacute;zier traits class is used
to represent: (i) an endpoint of a B&eacute;zier curve, (ii) a vertical tangency
point of a curve, used to subdivide it into \f$ x\f$-monotone subcurve, and
(iii) an intersection point between two curves. While, points of type (i) have
rational coordinates and are given as part of the input, points of the two
latter types have algebraic coordinates. However, to speed up the arrangement
construction, such point are not computed in an exact manner, and instead
are given in an approximate representation. Note that the exact coordinates
of a point may only be accessed if it is exactly computed.

In addition to the methods listed below, the copy constructor and assignment
operator for `Point_2` objects are also supported.

*/
class Point_2 {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Point_2 ();

/*!
constructs the point \f$ B(t_0)\f$ on the given curve. As \f$ t_0\f$ is an
algebraic number, the point has algebraic coordinates.
*/
Point_2 (const Curve_2& B, const Algebraic& t_0);

/*!
constructs the point \f$ B(t_0)\f$ on the given curve. As \f$ t_0\f$ is a
rational number, the point has rational coordinates.
*/
Point_2 (const Curve_2& B, const Rational& t_0);

/// @}

/// \name Access Functions
/// @{

/*!
returns the approximated coordinates of `p`.
*/
std::pair<double, double> approximate () const;

/*!
returns whether the coordinates of `p` are computed in an exact manner.
*/
bool is_exact () const;

/*!
returns the \f$ x\f$-coordinate of `p`.
\pre `p` is exactly computed.
*/
Algebraic x () const;

/*!
returns the \f$ y\f$-coordinate of `p`.
\pre `p` is exactly computed.
*/
Algebraic y () const;

/*!
returns whether the coordinates of `p` are rational numbers.
*/
bool is_rational () const;

/*!
casts `p` to a point with rational coordinates.
\pre `p` has rational coordinates.
*/
operator typename RatKernel::Point_2 () const;

/// @}

}; /* end Arr_Bezier_curve_traits_2::Point_2 */

/*!


The `X_monotone_curve_2` class nested within the B&eacute;zier traits is
used to represent \f$ x\f$-monotone subcurves of B&eacute;zier curves. The subcurve is
defined by a supporting B&eacute;zier curve \f$ B(t)\f$ and a range of definition in
the parameter space \f$ [t_1, t_2] \subseteq [0, 1]\f$, where \f$ B(t_1)\f$ is the
subcurve source and \f$ B(t_2)\f$ is its target. Note that as the point endpoints
may only be approximated, the parameter range defining the subcurve may
only be approximately known.

It is not possible to construct \f$ x\f$-monotone subcurves directly. Instead,
use the `Make_x_monotone_2` functor supplied by the traits class to
subdivide a `Curve_2` object into \f$ x\f$-monotone subcurves.

*/
class X_monotone_curve_2 {
public:

/// \name Access Functions
/// @{

/*!
returns the supporting B&eacute;zier curve of `b`.
*/
Curve_2 supporting_curve () const;

/*!
returns the source point of `b`.
*/
Point_2 source () const;

/*!
returns the target point of `b`.
*/
Point_2 target () const;

/*!
returns the left (\f$ xy\f$-lexicographically smaller) endpoint of `b`.
*/
Point_2 left () const;

/*!
returns the right (\f$ xy\f$-lexicographically smaller) endpoint of `b`.
*/
Point_2 right () const;

/*!
return the approximate parameter range defining the subcurve `b`.
*/
std::pair<double, double> parameter_range () const;

/// @}

}; /* end Arr_Bezier_curve_traits_2::X_monotone_curve_2 */

class Trim_2{
public:
/// \name Creation
/// @{

/*!
Trims the given x-monotone curve to an from src to tgt.
\ pre `src` and `tgt` lies on the curve
*/

X_monotone_curve_2(const X_monotone_curve_2& xcv,
                                const Point_2& src,
                                const Point_2& tgt)const

/// @}

}/* end Arr_Bezier_curve_traits_2::Trim_2 */

}; /* end Arr_Bezier_curve_traits_2 */
} /* end namespace CGAL */
