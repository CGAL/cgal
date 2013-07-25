
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` returns a const iterator range over all 
innermost coefficients of the given polynomial. 

\cgalRefines `AdaptableUnaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::ConstructInnermostCoefficientConstIteratorRange {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef std::pair< 
PolynomialTraits_d::Innermost_coefficient_const_iterator, 
PolynomialTraits_d::Innermost_coefficient_const_iterator > result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns a const iterator range over all innermost coefficients of \f$ p\f$. 
*/ 
result_type operator()(argument_type p); 

/// @}

}; /* end PolynomialTraits_d::ConstructInnermostCoefficientConstIteratorRange */

