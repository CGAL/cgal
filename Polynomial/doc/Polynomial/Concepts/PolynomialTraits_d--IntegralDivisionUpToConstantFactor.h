
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` computes the integral division 
of two polynomials of type `PolynomialTraits_d::Polynomial_d` 
<I>up to a constant factor (utcf)</I> . 

\pre \f$ g\f$ divides \f$ f\f$ in \f$ Q(R)[x_0,\dots,x_{d-1}]\f$, where \f$ Q(R)\f$ is the quotient field of the base ring \f$ R\f$, `PolynomialTraits_d::Innermost_coefficient_type`. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::GcdUpToConstantFactor`

*/

class PolynomialTraits_d::IntegralDivisionUpToConstantFactor {
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
Computes \f$ f/g\f$ up to a constant factor. 
*/ 
result_type operator()(first_argument_type f, 
second_argument_type g); 

/// @}

}; /* end PolynomialTraits_d::IntegralDivisionUpToConstantFactor */

