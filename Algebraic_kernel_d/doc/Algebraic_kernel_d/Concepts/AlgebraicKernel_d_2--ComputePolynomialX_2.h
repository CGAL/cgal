
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes a univariate square free polynomial \f$ p\f$, such that the first coordinate of
a given `AlgebraicKernel_d_2::Algebraic_real_2` is a real root of \f$ p\f$.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicKernel_d_2::ComputePolynomialY_2`

*/

class AlgebraicKernel_d_2::ComputePolynomialX_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_1 result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes a univariate square free polynomial \f$ p\f$, such that the first
coordinate of \f$ a\f$ is a real root of \f$ p\f$.
*/
result_type operator()(argument_type a);

/// @}

}; /* end AlgebraicKernel_d_2::ComputePolynomialX_2 */

