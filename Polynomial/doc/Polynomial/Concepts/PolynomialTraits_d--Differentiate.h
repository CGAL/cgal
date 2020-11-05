
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the derivative of a
`PolynomialTraits_d::Polynomial_d` with respect to one variable.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Differentiate {
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
Returns \f$ p'\f$, with respect to the outermost variable.
*/
result_type operator()(argument_type p);

/*!
Returns \f$ p'\f$, with respect to variable \f$ x_i\f$.
\pre \f$ 0 \leq i < d\f$.

*/
result_type operator()(argument_type p,
int i);

/// @}

}; /* end PolynomialTraits_d::Differentiate */

