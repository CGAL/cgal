
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes an open isolating interval for an `AlgebraicKernel_d_1::Algebraic_real_1`
with respect to the real roots of a given univariate polynomial.

\cgalRefines{AdaptableBinaryFunction}

\sa `AlgebraicKernel_d_1::ComputePolynomial_1`

*/

class AlgebraicKernel_d_1::Isolate_1 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::pair<AlgebraicKernel_d_1::Bound,AlgebraicKernel_d_1::Bound> result_type;

/*!

*/
typedef AlgebraicKernel_d_1::Algebraic_real_1 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_1::Polynomial_1 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes an open isolating interval \f$ I=(l,u)\f$ for \f$ a\f$ with respect to the real roots of \f$ p\f$.
It is not required that \f$ a\f$ is a root of \f$ p\f$.
\post \f$ a \in I\f$.
\post \f$ p(x) \neq0 | \forall x \in\overline{I}\backslash a\f$.

*/
result_type operator()(first_argument_type a, second_argument_type p);

/// @}

}; /* end AlgebraicKernel_d_1::Isolate_1 */

