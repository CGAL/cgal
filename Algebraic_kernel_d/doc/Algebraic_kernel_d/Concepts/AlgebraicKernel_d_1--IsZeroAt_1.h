
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Computes whether an `AlgebraicKernel_d_1::Polynomial_1` 
is zero at a given `AlgebraicKernel_d_1::Algebraic_real_1`. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_1::SignAt_1`

*/

class AlgebraicKernel_d_1::IsZeroAt_1 {
public:

/// \name Types 
/// @{

/*!
Type convertible to `bool` 
*/ 
typedef unspecified_type result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 first_argument_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes whether \f$ p\f$ is zero at \f$ x\f$. 
*/ 
result_type 
operator()(const first_argument_type & p, 
const second_argument_type & x); 

/// @}

}; /* end AlgebraicKernel_d_1::IsZeroAt_1 */

