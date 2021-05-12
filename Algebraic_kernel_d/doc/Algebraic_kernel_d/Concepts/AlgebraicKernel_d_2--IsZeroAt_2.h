
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes whether an `AlgebraicKernel_d_2::Polynomial_2`
is zero at a given `AlgebraicKernel_d_2::Algebraic_real_2`.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::SignAt_2`
\sa `AlgebraicKernel_d_1::IsZeroAt_1`

*/

class AlgebraicKernel_d_2::IsZeroAt_2 {
public:

/// \name Types
/// @{

/*!
Type convertible to `bool`
*/
typedef unspecified_type result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Polynomial_2 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes whether \f$ p\f$ is zero at \f$ a\f$.
*/
result_type
operator()(const first_argument_type & p,
const second_argument_type & a);

/// @}

}; /* end AlgebraicKernel_d_2::IsZeroAt_2 */

