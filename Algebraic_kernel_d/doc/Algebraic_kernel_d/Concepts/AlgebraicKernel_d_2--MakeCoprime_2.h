
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes for a given pair of bivariate polynomials \f$ p_1\f$, \f$ p_2\f$ their
common part \f$ g\f$ and coprime parts \f$ q_1\f$, \f$ q_2\f$ respectively.

That is, it computes \f$ g, q_1, q_2\f$ such that:

\f$ c_1 \cdot p_1 = g \cdot q_1\f$ for some constant \f$ c_1\f$ and

\f$ c_2 \cdot p_2 = g \cdot q_2\f$ for some constant \f$ c_2\f$,
such that \f$ q_1\f$ and \f$ q_2\f$ are coprime.

\cgalRefines `AdaptableFunctor` with five arguments

\sa `AlgebraicKernel_d_2::IsCoprime_2`

*/

class AlgebraicKernel_d_2::MakeCoprime_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef bool result_type;

/// @}

/// \name Operations
/// @{

/*!
Computes \f$ g, q_1, q_2\f$ as described above.

Returns whether \f$ p_1\f$ and \f$ p_2\f$ where already coprime.

*/
result_type
operator()(const AlgebraicKernel_d_2::Polynomial_2 & p1,
const AlgebraicKernel_d_2::Polynomial_2 & p2,
AlgebraicKernel_d_2::Polynomial_2 & g,
AlgebraicKernel_d_2::Polynomial_2 & q1,
AlgebraicKernel_d_2::Polynomial_2 & q2);

/// @}

}; /* end AlgebraicKernel_d_2::MakeCoprime_2 */

