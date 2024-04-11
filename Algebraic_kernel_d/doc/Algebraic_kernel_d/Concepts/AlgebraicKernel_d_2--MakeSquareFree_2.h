
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Returns a square free part of a bivariate polynomial.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicKernel_d_2::IsSquareFree_2`
\sa `AlgebraicKernel_d_2::SquareFreeFactorize_2`

*/

class AlgebraicKernel_d_2::MakeSquareFree_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 argument_type;

/// @}

/// \name Operations
/// @{

/*!

Returns a square free part of \f$ p\f$.
*/
result_type operator()(const argument_type& p);

/// @}

}; /* end AlgebraicKernel_d_2::MakeSquareFree_2 */

