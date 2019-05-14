
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Computes a square free univariate polynomial \f$ p\f$, such that the given 
`AlgebraicKernel_d_1::Algebraic_real_1` is a root of \f$ p\f$. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicKernel_d_1::Isolate_1`

*/

class AlgebraicKernel_d_1::ComputePolynomial_1 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes a square free polynomial \f$ p\f$, such that \f$ x\f$ is a real root of \f$ p\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end AlgebraicKernel_d_1::ComputePolynomial_1 */

