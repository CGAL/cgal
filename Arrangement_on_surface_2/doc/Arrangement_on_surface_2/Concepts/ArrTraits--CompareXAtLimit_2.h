namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableFunctor

\cgalHasModel ArrangementOpenBoundaryTraits_2::Compare_x_at_limit_2

*/

class CompareXAtLimit_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Given a point `p`, an \f$ x\f$-monotone curve `xcv`, and an
enumeration `ce` that specifies either the minimum or the
maximum end of the curve where the curve has a vertical asymptote,
compares the \f$ x\f$-coordinate of `p` and the \f$ x\f$-coordinate of the
limit of the curve at its specificed end. The variable `xcv`
identifies the parametric curve \f$ C(t) = (X(t),Y(t))\f$ defined over an
open or half-open interval with endpoints \f$ 0\f$ and \f$ 1\f$. The
enumeration `ce` identifies an open end \f$ d \in\{0,1\}\f$ of \f$ C\f$.
Formally, compares the \f$ x\f$-coordinate of `p` and
\f$ \lim_{t \rightarrow d} X(t)\f$. Returns `SMALLER`, `EQUAL`, or
`LARGER` accordingly.
\pre `parameter_space_in_y_2`(`xcv`, `ce`) \f$ \neq\f$ `ARR_INTERIOR`.
\pre If the parameter space is unbounded, \f$ C\f$ has a vertical asymptote at its \f$ d\f$-end; that is, `parameter_space_in_x_2`(`xcv`, `ce`) = `ARR_INTERIOR`.
*/
Comparison_result operator()(const ArrTraits::Point_2& p,
const ArrTraits::X_monotone_curve_2& xcv,
Arr_curve_end ce);

/*!
Given two \f$ x\f$-monotone curves `xcv1` and `xcv2` and two
indices `ce1` and `ce2` that specify either the minimum
or the maximum ends of `xcv1` and `xcv2`, respectively,
where the curves have vertical asymptotes, compares the
\f$ x\f$-coordinates of the limits of the curves at their specificed
ends. The variables `xcv1` and `xcv2` identify the
parametric curves \f$ C_1(t) = (X_1(t),Y_1(t))\f$ and
\f$ C_2(t) = (X_2(t),Y_2(t))\f$, respectively, defined over open or
half-open intervals with endpoints \f$ 0\f$ and \f$ 1\f$. The indices
`ce1` and `ce2` identify open ends \f$ d_1 \in\{0,1\}\f$ and
\f$ d_2 \in\{0,1\}\f$ of \f$ C_1\f$ and \f$ C_2\f$, respectively. Formally,
compares \f$ \lim_{t \rightarrow d_1} X_1(t)\f$ and
\f$ \lim_{t \rightarrow d_2} X_2(t)\f$. Returns `SMALLER`, `EQUAL`,
or `LARGER` accordingly.
\pre `parameter_space_in_y_2`(`xcv1`, `ce1`) \f$ \neq\f$ `ARR_INTERIOR`.
\pre `parameter_space_in_y_2`(`xcv2`, `ce2`) \f$ \neq\f$ `ARR_INTERIOR`.
\pre If the parameter space is unbounded, \f$ C_1\f$ has a vertical asymptote at its respective end; that is,
`parameter_space_in_x_2`(`xcv1`, `ce1`) = `ARR_INTERIOR`.
\pre If the parameter space is unbounded, \f$ C_2\f$ has a vertical asymptote at its respective end; that is,
`parameter_space_in_x_2`(`xcv2`, `ce2`) = `ARR_INTERIOR`.
*/
Comparison_result operator()(const ArrTraits::X_monotone_curve_2& xcv1,
Arr_curve_end ce1,
const ArrTraits::X_monotone_curve_2& xcv2,
Arr_curve_end ce2);

/// @}

}; /* end ArrTraits::CompareXAtLimit_2 */

}
