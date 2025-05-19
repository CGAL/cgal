
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the content of a
`PolynomialTraits_d::Polynomial_d`
with respect to the univariate (recursive) view on the
polynomial, that is, it computes the gcd of all
coefficients with respect to one variable.

This functor is well defined if `PolynomialTraits_d::Coefficient_type` is
a `Field` or a `UniqueFactorizationDomain`.

\cgalRefines{AdaptableUnaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::UnivariateContent {
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
Computes the content of \f$ p\f$ with respect to the outermost variable \f$ x_{d-1}\f$.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::UnivariateContent */

