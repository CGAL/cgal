
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes a square free factorization of an
`AlgebraicKernel_d_2::Polynomial_2`.

A polynomial \f$ p\f$ is factored into square free and pairwise
coprime non-constant factors \f$ q_i\f$ with multiplicities \f$ m_i\f$
and a constant factor \f$ c\f$, such that
\f$ p = c \cdot q_1^{m_1} \cdot ... \cdot q_n^{m_n}\f$.

The factor multiplicity pairs \f$ <q_i,m_i>\f$ are written to the
given output iterator. The constant factor \f$ c\f$ is not computed.

\cgalRefines{Assignable,CopyConstructible}

\sa `AlgebraicKernel_d_2::IsSquareFree_2`
\sa `AlgebraicKernel_d_2::MakeSquareFree_2`

*/

class AlgebraicKernel_d_2::SquareFreeFactorize_2 {
public:

/// \name Operations
/// @{

/*!
Copies in the output iterator the factors of a square free
factorization of \f$ p\f$, with their multiplicity, as objects of type
`std::pair<AlgebraicKernel_d_2::Polynomial_2, AlgebraicKernel_d_2::Multiplicity_type>`.
*/
template < class OutputIterator >
OutputIterator
operator()(const AlgebraicKernel_d_2::Polynomial_2& p, OutputIterator res);

/// @}

}; /* end AlgebraicKernel_d_2::SquareFreeFactorize_2 */

