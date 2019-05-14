
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Determines whether a given pair of univariate polynomials \f$ p_1, p_2\f$ is coprime, 
namely if \f$ \deg({\rm gcd}(p_1 ,p_2)) = 0\f$. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_1::MakeCoprime_1`

*/

class AlgebraicKernel_d_1::IsCoprime_1 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef bool result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 first_argument_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns true if `p1` and `p2` are coprime. 
*/ 
result_type 
operator()(const first_argument_type & p1, 
const second_argument_type & p2); 

/// @}

}; /* end AlgebraicKernel_d_1::IsCoprime_1 */

