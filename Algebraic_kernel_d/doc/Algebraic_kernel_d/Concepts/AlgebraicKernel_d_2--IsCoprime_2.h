
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes whether a given pair of bivariate polynomials is coprime.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::MakeCoprime_2`

*/

class AlgebraicKernel_d_2::IsCoprime_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef bool result_type;

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
Computes whether \f$ f\f$ and \f$ g\f$ are coprime.
*/
result_type
operator()(const first_argument_type& p1,
const second_argument_type& p2);

/// @}

}; /* end AlgebraicKernel_d_2::IsCoprime_2 */

