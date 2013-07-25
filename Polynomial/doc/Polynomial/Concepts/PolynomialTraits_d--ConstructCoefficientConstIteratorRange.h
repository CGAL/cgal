
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` returns a const iterator range over the 
coefficients of the given polynomial, with respect to the outermost variable, \f$ x_{d-1}\f$. 
The range starts with the coefficient for \f$ x_{d-1}^0\f$. 

\cgalRefines `AdaptableUnaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::ConstructCoefficientConstIteratorRange {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef std::pair< 
PolynomialTraits_d::Coefficient_const_iterator, 
PolynomialTraits_d::Coefficient_const_iterator > result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns a const iterator range over the coefficients of \f$ p\f$. 
*/ 
result_type operator()(argument_type p); 

/// @}

}; /* end PolynomialTraits_d::ConstructCoefficientConstIteratorRange */

