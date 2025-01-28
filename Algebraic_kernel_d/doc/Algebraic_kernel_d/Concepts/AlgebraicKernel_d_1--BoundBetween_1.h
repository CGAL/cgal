
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes a number of type
`AlgebraicKernel_d_1::Bound` in-between two
`AlgebraicKernel_d_1::Algebraic_real_1` values.

\cgalRefines{AdaptableBinaryFunction}

*/

class AlgebraicKernel_d_1::BoundBetween_1 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_1::Bound result_type;

/*!

*/
typedef AlgebraicKernel_d_1::Algebraic_real_1 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_1::Algebraic_real_1 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!

Computes a value \f$ r\f$, which is between \f$ a\f$ and \f$ b\f$.

\pre \f$ a \neq b\f$
\post \f$ r > min(a,b)\f$
\post \f$ r < max(a,b)\f$

*/
result_type
operator()(const first_argument_type & a,
const second_argument_type & b);

/// @}

}; /* end AlgebraicKernel_d_1::BoundBetween_1 */

