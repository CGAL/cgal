
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` computes the remainder of the 
<I>pseudo division</I> of two polynomials \f$ f\f$ and \f$ g\f$. 

Given \f$ f\f$ and \f$ g \neq 0\f$ one can compute quotient \f$ q\f$ and remainder \f$ r\f$ 
such that \f$ D \cdot f = g \cdot q + r\f$ and \f$ degree(r) < degree(g)\f$, 
where \f$ D = leading\_coefficient(g)^{max(0, degree(f)-degree(g)+1)}\f$ 

This functor computes \f$ r\f$. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::PseudoDivision`
\sa `PolynomialTraits_d::PseudoDivisionRemainder`
\sa `PolynomialTraits_d::PseudoDivisionQuotient`

*/

class PolynomialTraits_d::PseudoDivisionRemainder {
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

Returns the remainder \f$ r\f$ of the pseudo division of \f$ f\f$ and \f$ g\f$ with 
respect to the outermost variable \f$ x_{d-1}\f$. 
*/ 
result_type operator()(first_argument_type f, 
second_argument_type g); 

/// @}

}; /* end PolynomialTraits_d::PseudoDivisionRemainder */

