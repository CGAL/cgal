
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Computes the sign of a univariate polynomial 
`AlgebraicKernel_d_1::Polynomial_1` at a real value of type 
`AlgebraicKernel_d_1::Algebraic_real_1`. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_1::IsZeroAt_1`

*/

class AlgebraicKernel_d_1::SignAt_1 {
public:

/// \name Types 
/// @{

/*!
Type convertible to `CGAL::Sign` 
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
Computes the sign of \f$ p\f$ at \f$ x\f$. 
*/ 
result_type 
operator()(const first_argument_type & p, 
const second_argument_type & x); 

/// @}

}; /* end AlgebraicKernel_d_1::SignAt_1 */

