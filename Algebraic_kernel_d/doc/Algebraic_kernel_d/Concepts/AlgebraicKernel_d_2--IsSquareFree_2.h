
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes whether the given bivariate polynomial is square free.

\cgalRefines `AdaptableUnaryFunction`

\sa `AlgebraicKernel_d_2::MakeSquareFree_2`
\sa `AlgebraicKernel_d_2::SquareFreeFactorize_2`

*/

class AlgebraicKernel_d_2::IsSquareFree_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef bool result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 argument_type;

/// @}

/// \name Operations
/// @{

/*!

Computes whether \f$ p\f$ is square free.
*/
result_type operator()(const argument_type& p);

/// @}

}; /* end AlgebraicKernel_d_2::IsSquareFree_2 */

