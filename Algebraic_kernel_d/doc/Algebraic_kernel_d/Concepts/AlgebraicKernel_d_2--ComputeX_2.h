
/*!
\ingroup PkgAlgebraicKerneldConceptsBi
\cgalConcept

Computes the first coordinate of an 
`AlgebraicKernel_d_2::AlgebraicReal_2`. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicKernel_d_2::ComputeY_2`

*/

class AlgebraicKernel_d_2::ComputeX_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_2::Algebraic_real_1 result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_2::Algebraic_real_2 argument_type; 

/// @} 

/// \name Operations 
/// A model of this type must provide:
/// @{

/*!

Computes the first coordinate of \f$ a\f$. 

*/ 
result_type operator()(const argument_type & a); 

/// @}

}; /* end AlgebraicKernel_d_2::ComputeX_2 */

