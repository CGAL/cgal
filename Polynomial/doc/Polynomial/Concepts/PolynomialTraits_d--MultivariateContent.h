
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the content of a
`PolynomialTraits_d::Polynomial_d` with respect to the symmetric
view on the polynomial, that is, it computes the gcd of all innermost coefficients.

This functor is well defined if `PolynomialTraits_d::Innermost_coefficient_type` is a
`Field` or a `UniqueFactorizationDomain`.

\cgalRefines{AdaptableUnaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::MultivariateContent {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Innermost_coefficient_type result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes the \f$ gcd\f$ of all innermost coefficients of \f$ p\f$.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::MultivariateContent */

