
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes the real zero-dimensional solutions of a bivariate polynomial system.
The multiplicity stored in the output iterator is the multiplicity in the system.

\cgalRefines `Assignable`
\cgalRefines `CopyConstructible`
*/
class AlgebraicKernel_d_2::Solve_2 {
public:

/// \name Operations
/// A model of this type must provide:
/// @{

/*!
Computes all common solutions of \f$ f\f$ and \f$ g\f$ with multiplicity, and copies them as objects of type
`std::pair<AlgebraicKernel_d_2::Algebraic_real_2, AlgebraicKernel_d_2::Multiplicity_type>` in `res`.
\pre \f$ f\f$ is square free.
\pre \f$ g\f$ is square free.
\pre \f$ f\f$ and \f$ g\f$ are coprime.

*/
template < class OutputIterator >
OutputIterator
operator()(
AlgebraicKernel_d_2::Polynomial_2 f,
AlgebraicKernel_d_2::Polynomial_2 g,
OutputIterator res );

/*!
Computes all common solutions of \f$ f\f$ and \f$ g\f$ in the closed box \f$ [xl,xu]\times[yl,yu]\f$, and copies them as objects of
type `std::pair<AlgebraicKernel_d_2::Algebraic_real_2, AlgebraicKernel_d_2::Multiplicity_type>` in `res`.
\pre \f$ f\f$ is square free.
\pre \f$ g\f$ is square free.
\pre \f$ f\f$ and \f$ g\f$ are coprime.

*/
template < class OutputIterator > OutputIterator operator()(
AlgebraicKernel_d_2::Polynomial_2 f,
AlgebraicKernel_d_2::Polynomial_2 g,
AlgebraicKernel_d_2::Bound xl,
AlgebraicKernel_d_2::Bound xu,
AlgebraicKernel_d_2::Bound yl,
AlgebraicKernel_d_2::Bound yu,
OutputIterator res);

/// @}

}; /* end AlgebraicKernel_d_2::Solve_2 */

