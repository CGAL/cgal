
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes an isolating interval for the first coordinate of an `AlgebraicKernel_d_2::Algebraic_real_2`
with respect to the real roots of a univariate polynomial.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::IsolateY_2`
\sa `AlgebraicKernel_d_2::ComputePolynomialX_2`
\sa `AlgebraicKernel_d_2::ComputePolynomialY_2`

*/

class AlgebraicKernel_d_2::IsolateX_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::pair<AlgebraicKernel_d_2::Bound,AlgebraicKernel_d_2::Bound> result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_1 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes an open isolating interval \f$ I=(l,u)\f$ for the first coordinate \f$ x\f$ of \f$ a\f$ with respect to the real roots of \f$ p\f$.
It is not required that \f$ x\f$ is a root of \f$ p\f$.
\post \f$ x \in I\f$.
\post \f$ p(\alpha) \neq0 | \forall\alpha\in\overline{I}\backslash x\f$.

*/
result_type operator()(first_argument_type a, second_argument_type p);

/// @}

}; /* end AlgebraicKernel_d_2::IsolateX_2 */

