
/*!
\ingroup PkgPolynomialConcepts
\cgalconcept

This `AdaptableUnaryFunction` computes whether the given 
a polynomial of type `PolynomialTraits_d::Polynomial_d` 
is square free. 

Note that this statement does cover constant factors, 
i.e., whether the multivariate content contains a square. 

\refines ::AdaptableUnaryFunction 
\refines ::CopyConstructible 
\refines ::DefaultConstructible 

\sa  \ref ::Polynomial_d 
\sa  \ref ::PolynomialTraits_d 
\sa  \ref ::PolynomialTraits_d::SquareFreeFactorize 
\sa  \ref ::PolynomialTraits_d::MakeSquareFree 
\sa  \ref ::PolynomialTraits_d::MultivariateContent 

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

