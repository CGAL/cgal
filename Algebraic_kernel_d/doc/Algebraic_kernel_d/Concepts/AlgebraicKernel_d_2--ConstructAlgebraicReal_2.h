
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Constructs an `AlgebraicKernel_d_2::Algebraic_real_2`.

\cgalRefines `AdaptableFunctor`

\sa `AlgebraicKernel_d_1::ConstructAlgebraicReal_1`

*/

class AlgebraicKernel_d_2::ConstructAlgebraicReal_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 result_type;

/// @}

/// \name Operations
/// @{

/*!
introduces an `AlgebraicKernel_d_2::Algebraic_real_2` initialized to \f$ (x,y)\f$.
*/
result_type operator()(int x, int y);

/*!
introduces an `AlgebraicKernel_d_2::Algebraic_real_2` initialized to \f$ (x,y)\f$.
*/
result_type operator()(AlgebraicKernel_d_2::Bound x,AlgebraicKernel_d_2::Bound y);

/*!
introduces an `AlgebraicKernel_d_2::Algebraic_real_2` initialized to \f$ (x,y)\f$.
*/
result_type operator()(AlgebraicKernel_d_2::Coefficient x,AlgebraicKernel_d_2::Coefficient y);

/*!
introduces an `AlgebraicKernel_d_2::Algebraic_real_2` initialized to \f$ (x,y)\f$.
*/
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_1 x,AlgebraicKernel_d_2::Algebraic_real_1 y);

/*!

introduces an `AlgebraicKernel_d_2::Algebraic_real_2`
initialized to the \f$ i\f$-th real common solution of \f$ f\f$ and \f$ g\f$,
with respect to `xy`-lexicographic order.
The index starts at \f$ 0\f$, that is, the system must have at least \f$ i+1\f$ real solutions.
\pre \f$ f\f$ is square free.
\pre \f$ g\f$ is square free.
\pre \f$ f\f$ and \f$ g\f$ are coprime.

*/
result_type operator()(
AlgebraicKernel_d_2::Polynomial_2 f,
AlgebraicKernel_d_2::Polynomial_2 g,
AlgebraicKernel_d_2::size_type i);

/*!
introduces an `AlgebraicKernel_d_2::Algebraic_real_2` initialized to the only real intersection of
\f$ f\f$ and \f$ g\f$ in the open box \f$ B = (x_l,x_u)\times(y_l,y_u)\f$.
\pre \f$ x_l < x_u\f$
\pre \f$ y_l < y_u\f$
\pre \f$ f\f$ is square free.
\pre \f$ g\f$ is square free.
\pre \f$ f\f$ and \f$ g\f$ are coprime.
\pre \f$ f\f$ and \f$ g\f$ have exactly one common solution in \f$ B\f$
\pre \f$ f\f$ and \f$ g\f$ have no common solution on \f$ \partial B\f$

*/
result_type operator()(
AlgebraicKernel_d_2::Polynomial_2 f,
AlgebraicKernel_d_2::Polynomial_2 g,
AlgebraicKernel_d_2::Bound x_l,
AlgebraicKernel_d_2::Bound x_u,
AlgebraicKernel_d_2::Bound y_l,
AlgebraicKernel_d_2::Bound y_u);

/// @}

}; /* end AlgebraicKernel_d_2::ConstructAlgebraicReal_2 */

