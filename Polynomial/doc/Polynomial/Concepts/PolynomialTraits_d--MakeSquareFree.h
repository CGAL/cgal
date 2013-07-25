
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the square-free part of 
a polynomial of type `PolynomialTraits_d::Polynomial_d` 
<I>up to a constant factor</I>. 

A polynomial \f$ p\f$ can be factored into square-free and pairwise coprime 
non-constant factors \f$ g_i\f$ with multiplicities \f$ m_i\f$ and a constant factor \f$ a\f$, 
such that \f$ p = a \cdot g_1^{m_1} \cdot ... \cdot g_n^{m_n}\f$, where all \f$ g_i\f$ are canonicalized. 

Given this decomposition, the square free part is defined as the product \f$ g_1 \cdot ... \cdot g_n\f$, 
which is computed by this functor. 

\cgalRefines `AdaptableUnaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Canonicalize`
\sa `PolynomialTraits_d::SquareFreeFactorize`
\sa `PolynomialTraits_d::IsSquareFree`

*/

class PolynomialTraits_d::MakeSquareFree {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the square-free part of \f$ p\f$. 
*/ 
result_type operator()(argument_type p); 

/// @}

}; /* end PolynomialTraits_d::MakeSquareFree */

