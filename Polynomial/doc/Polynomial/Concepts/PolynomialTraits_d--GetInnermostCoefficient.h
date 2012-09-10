
/*!
\ingroup PkgPolynomialConcepts
\cgalconcept

For the given `PolynomialTraits_d::Polynomial_d` this 
`AdaptableBinaryFunction` returns the coefficient of 
the (multivariate) monomial specified by the given `Exponent_vector`. 

\refines ::AdaptableBinaryFunction 
\refines ::CopyConstructible 
\refines ::DefaultConstructible 

\sa  \ref ::Polynomial_d 
\sa  \ref ::PolynomialTraits_d 

*/

class PolynomialTraits_d::GetInnermostCoefficient {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef PolynomialTraits_d::Innermost_coefficient_type result_type; 

/*! 

*/ 
typedef PolynomialTraits_d::Polynomial_d first_argument_type ; 

/*! 

*/ 
typedef Exponent_vector second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 

For given polynomial \f$ p\f$ this operator returns the innermost coefficient of the 
monomial corresponding to the given `Exponent_vector` \f$ v\f$. 
*/ 
result_type operator()( first_argument_type p, 
second_argument_type v); 

/// @}

}; /* end PolynomialTraits_d::GetInnermostCoefficient */

