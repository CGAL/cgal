
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes an isolating box for a given `AlgebraicKernel_d_2::Algebraic_real_2`.

\cgalRefines{AdaptableFunctor}

\sa `AlgebraicKernel_d_2::IsolateX_2`
\sa `AlgebraicKernel_d_2::IsolateY_2`
\sa `AlgebraicKernel_d_2::ComputePolynomialX_2`
\sa `AlgebraicKernel_d_2::ComputePolynomialY_2`

*/

class AlgebraicKernel_d_2::Isolate_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::array<AlgebraicKernel_d_1::Bound, 4> result_type;

/// @}

/// \name Operations
/// @{

/*!
The returned `std::array` \f$ [xl,xu,yl,yu]\f$ represents an open isolating box \f$ B=(xl,xu)\times(yl,yu)\f$
for \f$ a\f$ with respect to \f$ f\f$.
\pre \f$ f(a)\neq0\f$
\post \f$ a \in B\f$.
\post \f$ \{ r | f(r)=0 \} \cap\overline{B} = \emptyset\f$.

*/
result_type
operator()( AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Polynomial_2 f);

/*!
The returned `std::array` \f$ [xl,xu,yl,yu]\f$ represents an open isolating box \f$ B=(xl,xu)\times(yl,yu)\f$
for \f$ a\f$ with respect to the common solutions of \f$ f\f$ and \f$ g\f$.
It is not necessary that \f$ a\f$ is a common solution of \f$ f\f$ and \f$ g\f$.
\post \f$ a \in B\f$.
\post \f$ \{ r | f(r)=g(r)=0 \} \cap\overline{B} \in\{\{a\},\emptyset\}\f$.

*/
result_type
operator()(
AlgebraicKernel_d_2::Algebraic_real_2 a,
AlgebraicKernel_d_2::Polynomial_2 f,
AlgebraicKernel_d_2::Polynomial_2 g);

/// @}

}; /* end AlgebraicKernel_d_2::Isolate_2 */

