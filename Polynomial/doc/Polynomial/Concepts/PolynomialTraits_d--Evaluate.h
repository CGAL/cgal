
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` evaluates 
`PolynomialTraits_d::Polynomial_d` with respect to one variable. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Evaluate {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef PolynomialTraits_d::Coefficient_type result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d first_argument_type; 

/*!

*/ 
typedef PolynomialTraits_d::Coefficient_type second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns \f$ p(x)\f$, with respect to the outermost variable. 
*/ 
result_type operator()(first_argument_type p, 
second_argument_type x); 

/// @}

}; /* end PolynomialTraits_d::Evaluate */

