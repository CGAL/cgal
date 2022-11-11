
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes whether the given
a polynomial of type `PolynomialTraits_d::Polynomial_d`
is square free.

Note that this statement does cover constant factors,
i.e., whether the multivariate content contains a square.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::SquareFreeFactorize`
\sa `PolynomialTraits_d::MakeSquareFree`
\sa `PolynomialTraits_d::MultivariateContent`

*/

class PolynomialTraits_d::IsSquareFree {
public:

/// \name Types
/// @{

/*!

*/
typedef bool result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns whether the \f$ p\f$ is square free.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::IsSquareFree */

