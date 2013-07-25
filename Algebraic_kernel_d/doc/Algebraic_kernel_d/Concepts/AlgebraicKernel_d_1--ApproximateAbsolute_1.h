
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

A model of `AlgebraicKernel_d_1::ApproximateAbsolute_1` is an `AdaptableBinaryFunction` that computes an 
approximation of an `AlgebraicKernel_d_1::Algebraic_real_1` value with 
respect to a given absolute precision. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_1::ApproximateRelative_1`

*/

class AlgebraicKernel_d_1::ApproximateAbsolute_1 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef std::pair<AlgebraicKernel_d_1::Bound, AlgebraicKernel_d_1::Bound> result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 first_argument_type; 

/*!

*/ 
typedef int second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!

The function computes a pair \f$ p\f$ of `AlgebraicKernel_d_1::Bound`, 
where \f$ p.first\f$ represents the lower approximation and \f$ p.second\f$ represents 
the upper approximation. The pair \f$ p\f$ approximates the given value \f$ x\f$ 
with respect to the given absolute precision \f$ a\f$. 
\post \f$ p.first <= x \f$ 
\post \f$ x <= p.second \f$ 
\post \f$ (x - p.first) <= 2^{-a} \f$ 
\post \f$ (p.second - x) <= 2^{-a} \f$ 

*/ 
result_type 
operator()(const first_argument_type & x, 
const second_argument_type & a ); 

/// @}

}; /* end AlgebraicKernel_d_1::ApproximateAbsolute_1 */

