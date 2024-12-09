
/*!
\ingroup PkgAlgebraicKernelDConceptsUni
\cgalConcept

Computes the number of real solutions of the given univariate polynomial.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicKernel_d_1::ConstructAlgebraicReal_1`

*/

class AlgebraicKernel_d_1::NumberOfSolutions_1 {
public:

/// \name Types
/// A model of this type must provide:
/// @{

/*!

*/
typedef AlgebraicKernel_d_1::size_type result_type;

/*!

*/
typedef AlgebraicKernel_d_1::Polynomial_1 argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns the number of real solutions of \f$ p\f$.
\pre \f$ p\f$ is square free.
*/
result_type operator()(argument_type p);

/// @}

}; /* end AlgebraicKernel_d_1::NumberOfSolutions_1 */

