
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` translates a
`PolynomialTraits_d::Polynomial_d` with respect to one variable, that is,
for a given polynomial \f$ p\f$ and constant \f$ c\f$ it computes \f$ p(x+c)\f$.

Note that this functor operates on the polynomial in the univariate view, that is,
the polynomial is considered as a univariate polynomial in one specific variable.

\cgalRefines `AdaptableBinaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Translate {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Polynomial_d result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d first_argument_type;

/*!

*/
typedef PolynomialTraits_d::Innermost_coefficient_type second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns \f$ p(x+c)\f$, with respect to the outermost variable.
*/
result_type operator()(first_argument_type p,
second_argument_type c);

/*!
Same as first operator but for variable \f$ x_i\f$.
\pre \f$ 0 \leq i < d\f$.

*/
result_type operator()(first_argument_type p,
second_argument_type c,
int i);

/// @}

}; /* end PolynomialTraits_d::Translate */

