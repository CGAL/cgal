
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the innermost leading coefficient
of a `PolynomialTraits_d::Polynomial_d` \f$ p\f$. The innermost leading coefficient is recursively defined as the innermost leading coefficient of the leading coefficient of \f$ p\f$. In case \f$ p\f$ is univariate it coincides with the leading coefficient.

\cgalRefines{AdaptableUnaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::InnermostLeadingCoefficient {
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
Computes the innermost leading coefficient of \f$ p\f$.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::InnermostLeadingCoefficient */

