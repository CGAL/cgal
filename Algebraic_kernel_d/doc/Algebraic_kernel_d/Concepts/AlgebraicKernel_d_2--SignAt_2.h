
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes the sign of a bivariate polynomial
`AlgebraicKernel_d_2::Polynomial_2` at a value of type
`AlgebraicKernel_d_2::Algebraic_real_2`.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::IsZeroAt_2`
\sa `AlgebraicKernel_d_1::SignAt_1`

*/

class AlgebraicKernel_d_2::SignAt_2 {
public:

/// \name Types
/// @{

/*!
Type convertible to `CGAL::Sign`
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
Computes the sign of a bivariate polynomial \f$ p\f$ evaluated at \f$ a\f$.
*/
result_type
operator()(const first_argument_type & p, const second_argument_type & a);

/// @}

}; /* end AlgebraicKernel_d_2::SignAt_2 */

