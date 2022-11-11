
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the content of a
`PolynomialTraits_d::Polynomial_d`
with respect to the univariate (recursive) view on the
polynomial <I>up to a constant factor (utcf)</I>, that is,
it computes the \f$ \mathrm{gcd\_utcf}\f$ of all coefficients with respect to one variable.

Remark: This is called `UnivariateContentUpToConstantFactor` for
symmetric reasons with respect to `PolynomialTraits_d::UnivariateContent`
and `PolynomialTraits_d::MultivariateContent`.
However, a concept `PolynomialTraits_d::MultivariateContentUpToConstantFactor`
does not exist since the result is trivial.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::GcdUpToConstantFactor`

*/

class PolynomialTraits_d::UnivariateContentUpToConstantFactor {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Coefficient_type result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes the content <I>up to a constant factor</I> of \f$ p\f$ with
respect to the outermost variable \f$ x_{d-1}\f$.
*/
result_type operator()(first_argument_type p);

/// @}

}; /* end PolynomialTraits_d::UnivariateContentUpToConstantFactor */

