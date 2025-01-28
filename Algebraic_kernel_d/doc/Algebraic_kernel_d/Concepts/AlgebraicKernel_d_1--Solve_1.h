
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes the real roots of a univariate polynomial.

\cgalRefines{Assignable,CopyConstructible}

*/

class AlgebraicKernel_d_1::Solve_1 {
public:

/// \name Operations
/// A model of this type must provide:
/// @{

/*!
Computes all real solutions of \f$ p\f$ with multiplicity, and copies them as objects of type
`std::pair<AlgebraicKernel_d_1::Algebraic_real_1, AlgebraicKernel_d_1::Multiplicity_type>` in `res`.
*/
template < class OutputIterator >
OutputIterator
operator()( AlgebraicKernel_d_1::Polynomial_1 p,
OutputIterator res);

/*!
Computes all real solutions of \f$ p\f$, and copies them as objects of type
`AlgebraicKernel_d_1::Algebraic_real_1` in `res`.
The `bool` `known_to_be_square_free` indicates whether \f$ p\f$ is known to be square free.
Each root, though it might be a multiple root, is reported only once.
*/
template < class OutputIterator >
OutputIterator
operator()(
AlgebraicKernel_d_1::Polynomial_1 p,
bool known_to_be_square_free,
OutputIterator res );

/*!
Computes all real solutions of \f$ p\f$ in the closed interval \f$ [l,u]\f$ with multiplicity, and copies them as objects of type
`std::pair<AlgebraicKernel_d_1::Algebraic_real_1, AlgebraicKernel_d_1::Multiplicity_type>` in `res`.
*/
template < class OutputIterator >
OutputIterator
operator()( AlgebraicKernel_d_1::Polynomial_1 p, AlgebraicKernel_d_1::Bound l, AlgebraicKernel_d_1::Bound u,
OutputIterator res);

/*!
Computes all real solutions of \f$ p\f$ in the closed interval \f$ [l,u]\f$, and copies them as objects of
type `AlgebraicKernel_d_1::Algebraic_real_1` in `res`.
The `bool` `known_to_be_square_free` indicates whether \f$ p\f$ is known to be square free.
Each root, though it might be a multiple root, is reported only once.
*/
template < class OutputIterator > OutputIterator operator()(
AlgebraicKernel_d_1::Polynomial_1 p,
bool known_to_be_square_free,
AlgebraicKernel_d_1::Bound l,
AlgebraicKernel_d_1::Bound u,
OutputIterator res);

/// @}

}; /* end AlgebraicKernel_d_1::Solve_1 */

