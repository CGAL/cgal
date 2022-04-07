
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes the number of real solutions of the given bivariate polynomial system.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::ConstructAlgebraicReal_2`

*/

class AlgebraicKernel_d_2::NumberOfSolutions_2 {
public:

/// \name Types
/// A model of this type must provide:
/// @{

/*!

*/
typedef AlgebraicKernel_d_2::size_type result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns the number of real solutions of the bivariate polynomial system \f$ (f,g)\f$.
\pre \f$ f\f$ is square free.
\pre \f$ g\f$ is square free.
\pre \f$ f\f$ and \f$ g\f$ are coprime.
*/
result_type operator()(first_argument_type f, second_argument_type g);

/// @}

}; /* end AlgebraicKernel_d_2::NumberOfSolutions_2 */

