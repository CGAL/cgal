
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes \f$ p(-x)\f$ for a given polynomial \f$ p\f$.

Note that this functor operates on the polynomial in the univariate view, that is,
the polynomial is considered as a univariate polynomial in one specific variable.

This functor is provided for efficiency reasons, since this operation just flips the sign
of all odd coefficients with respect to the specified variable.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Negate {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Polynomial_d result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns \f$ p(-x)\f$, with respect to the outermost variable.
*/
result_type operator()(argument_type p);

/*!
Returns \f$ p(-x)\f$, with respect to variable \f$ x_i\f$.
\pre \f$ 0 \leq i < d\f$.

*/
result_type operator()(argument_type p,
int i);

/// @}

}; /* end PolynomialTraits_d::Negate */

