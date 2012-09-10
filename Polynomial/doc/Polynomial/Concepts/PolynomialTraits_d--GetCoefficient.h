
/*!
\ingroup PkgPolynomialConcepts
\cgalconcept

This `AdaptableBinaryFunction` provides access to coefficients of a 
`PolynomialTraits_d::Polynomial_d`. 

\refines ::AdaptableBinaryFunction 
\refines ::CopyConstructible 
\refines ::DefaultConstructible 

\sa  \ref ::Polynomial_d 
\sa  \ref ::PolynomialTraits_d 

*/

class PolynomialTraits_d::GetCoefficient {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef PolynomialTraits_d::Coefficient_type result_type; 

/*! 

*/ 
typedef PolynomialTraits_d::Polynomial_d first_argument_type ; 

/*! 

*/ 
typedef int second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 

For given polynomial \f$ p\f$ this operator returns the coefficient 
of \f$ x_{d-1}^e\f$, where \f$ x_{d-1}\f$ is the outermost variable. 
*/ 
result_type operator()( first_argument_type p, 
second_argument_type e); 

/// @}

}; /* end PolynomialTraits_d::GetCoefficient */

