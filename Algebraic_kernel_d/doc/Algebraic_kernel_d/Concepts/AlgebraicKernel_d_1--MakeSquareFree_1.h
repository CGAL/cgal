
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Returns a square free part of a univariate polynomial.

\cgalRefines `AdaptableUnaryFunction`

\sa `AlgebraicKernel_d_1::IsSquareFree_1`
\sa `AlgebraicKernel_d_1::SquareFreeFactorize_1`

*/

class AlgebraicKernel_d_1::MakeSquareFree_1 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_1::Polynomial_1 result_type;

/*!

*/
typedef AlgebraicKernel_d_1::Polynomial_1 argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns a square free part of \f$ p\f$
*/
result_type operator()(argument_type p);

/// @}

}; /* end AlgebraicKernel_d_1::MakeSquareFree_1 */

