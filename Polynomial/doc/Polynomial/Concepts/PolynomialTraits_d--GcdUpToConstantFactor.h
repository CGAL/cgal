
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` computes the \f$ gcd\f$ 
<I>up to a constant factor (utcf)</I> of two polynomials of type 
`PolynomialTraits_d::Polynomial_d`. 

In case the base ring \f$ R\f$ (`PolynomialTraits_d::Innermost_coefficient_type`) 
is not a `UniqueFactorizationDomain` or not a `Field` the polynomial ring 
\f$ R[x_0,\dots,x_{d-1}]\f$ (`PolynomialTraits_d::Polynomial_d`) may not 
possesses greatest common divisors. However, since \f$ R\f$ is an integral 
domain one can consider its quotient field \f$ Q(R)\f$ for which \f$ gcd\f$s of 
polynomials exist. 

This functor computes \f$ gcd\_utcf(f,g) = D * gcd(f,g)\f$, 
for some \f$ D \in R\f$ such that \f$ gcd\_utcf(f,g) \in R[x_0,\dots,x_{d-1}]\f$. 
Hence, \f$ gcd\_utcf(f,g)\f$ may not be a divisor of \f$ f\f$ and \f$ g\f$ in \f$ R[x_0,\dots,x_{d-1}]\f$. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::IntegralDivisionUpToConstantFactor`
\sa `PolynomialTraits_d::UnivariateContentUpToConstantFactor`
\sa `PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor`

*/

class PolynomialTraits_d::GcdUpToConstantFactor {
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
typedef PolynomialTraits_d::Polynomial_d second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes \f$ gcd(f,g)\f$ up to a constant factor. 
*/ 
result_type operator()(first_argument_type f, 
second_argument_type g); 

/// @}

}; /* end PolynomialTraits_d::GcdUpToConstantFactor */

