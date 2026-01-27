
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes for a given pair of univariate polynomials \f$ p_1\f$, \f$ p_2\f$ their
common part \f$ g\f$ up to a constant factor and coprime parts \f$ q_1\f$, \f$ q_2\f$
respectively.

That is, it computes \f$ g, q_1, q_2\f$ such that:

\f$ c_1 \cdot p_1 = g \cdot q_1\f$ for some constant \f$ c_1\f$ and

\f$ c_2 \cdot p_2 = g \cdot q_2\f$ for some constant \f$ c_2\f$,
such that \f$ q_1\f$ and \f$ q_2\f$ are coprime.

It returns true if \f$ p_1\f$ and \f$ p_2\f$ are already coprime.

\cgalRefines{AdaptableQuinaryFunction}

\sa `AlgebraicKernel_d_1::IsCoprime_1`

*/

class AlgebraicKernel_d_1::MakeCoprime_1 {
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
operator()(const AlgebraicKernel_d_1::Polynomial_1 & p1,
const AlgebraicKernel_d_1::Polynomial_1 & p2,
AlgebraicKernel_d_1::Polynomial_1 & g,
AlgebraicKernel_d_1::Polynomial_1 & q1,
AlgebraicKernel_d_1::Polynomial_1 & q2);

/// @}

}; /* end AlgebraicKernel_d_1::MakeCoprime_1 */

