
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` multiplies a `PolynomialTraits_d::Polynomial_d` 
by the given power of the specified variable. 

This functor is provided for efficiency reasons, since multiplication by some variable 
will in general correspond to a shift of coefficients in the internal representation. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Shift {
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
typedef int second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns \f$ p * x_{d-1}^e\f$. 
\pre \f$ 0 \leq e\f$. 
*/ 
result_type operator()(first_argument_type p, 
second_argument_type e); 

/*!
Returns \f$ p * x_{i}^e\f$. 
\pre \f$ 0 \leq e\f$. 
\pre \f$ 0 \leq i < d\f$. 

*/ 
result_type operator()(first_argument_type p, 
second_argument_type e, 
int i); 

/// @}

}; /* end PolynomialTraits_d::Shift */

