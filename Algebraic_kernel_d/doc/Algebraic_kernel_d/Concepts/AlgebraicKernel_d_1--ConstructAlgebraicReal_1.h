
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Constructs `AlgebraicKernel_d_1::Algebraic_real_1`. 

\cgalRefines `AdaptableFunctor` 

\sa `AlgebraicKernel_d_2::ConstructAlgebraicReal_2`

*/

class AlgebraicKernel_d_1::ConstructAlgebraicReal_1 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 result_type; 

/// @} 

/// \name Operations 
/// @{

/*!
introduces an `AlgebraicKernel_d_1::Algebraic_real_1` initialized to \f$ a\f$. 
*/ 
result_type operator()(int a); 

/*!
introduces an `AlgebraicKernel_d_1::Algebraic_real_1` initialized to \f$ a\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Bound a); 

/*!
introduces an `AlgebraicKernel_d_1::Algebraic_real_1` initialized to \f$ a\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Coefficient a); 

/*!
introduces an `AlgebraicKernel_d_1::Algebraic_real_1` initialized to the \f$ i\f$-th real root of \f$ p\f$. 
The index starts at \f$ 0\f$, that is, \f$ p\f$ must have at least \f$ i+1\f$ real roots. 
\pre \f$ p\f$ is square free. 
\pre \f$ p\f$ has at least \f$ i+1\f$ real roots. 

*/ 
result_type operator()(AlgebraicKernel_d_1::Polynomial_1 p, AlgebraicKernel_d_1::size_type i); 

/*!
introduces an `AlgebraicKernel_d_1::Algebraic_real_1` initialized to the 
only real root of \f$ p\f$ in the open interval \f$ I = (l,u)\f$. 
\pre \f$ l < u\f$ 
\pre \f$ p\f$ is square free. 
\pre \f$ p\f$ has exactly one real root in \f$ I\f$ 
\pre \f$ p\f$ has no real root on \f$ \partial I\f$ 

*/ 
result_type operator()( 
AlgebraicKernel_d_1::Polynomial_1 p, 
AlgebraicKernel_d_1::Bound l, 
AlgebraicKernel_d_1::Bound u); 

/// @}

}; /* end AlgebraicKernel_d_1::ConstructAlgebraicReal_1 */

