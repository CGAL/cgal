
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the leading coefficient 
of a `PolynomialTraits_d::Polynomial_d`. 

\cgalRefines `AdaptableUnaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::LeadingCoefficient {
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

Computes the leading coefficient of \f$ p\f$ with respect to the 
outermost variable \f$ x_{d-1}\f$. 
*/ 
result_type operator()(argument_type p); 

/// @}

}; /* end PolynomialTraits_d::LeadingCoefficient */

