
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes whether the given univariate polynomial is square free.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicKernel_d_1::MakeSquareFree_1`
\sa `AlgebraicKernel_d_1::SquareFreeFactorize_1`

*/

class AlgebraicKernel_d_1::IsSquareFree_1 {
public:

/// \name Types
/// A model of this type must provide:
/// @{

/*!

*/
typedef bool result_type;

/*!

*/
typedef AlgebraicKernel_d_1::Polynomial_1 argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns true if \f$ p\f$ is square free.
*/
result_type operator()(argument_type p);

/// @}

}; /* end AlgebraicKernel_d_1::IsSquareFree_1 */

