
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Compares `AlgebraicKernel_d_2::Algebraic_real_2`s lexicographically.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicKernel_d_2::CompareX_2`
\sa `AlgebraicKernel_d_2::CompareY_2`

*/

class AlgebraicKernel_d_2::CompareXY_2 {
public:

/// \name Types
/// @{

/*!
Type convertible to `CGAL::Comparison_result`
*/
typedef unspecified_type result_type;

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 first_argument_type;

/*!

*/
typedef AlgebraicKernel_d_2::Algebraic_real_2 second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Compares \f$ a\f$ and \f$ b\f$ lexicographically.
*/
result_type operator()(const first_argument_type & a, const second_argument_type & b);

/*!
Compares \f$ a\f$ with \f$ (x,y)\f$ lexicographically.
*/
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, int x, int y);

/*!
Compares \f$ a\f$ with \f$ (x,y)\f$ lexicographically.
*/
result_type operator()(
AlgebraicKernel_d_2::Algebraic_real_2 a,
AlgebraicKernel_d_2::Bound x,
AlgebraicKernel_d_2::Bound y);

/*!
Compares \f$ a\f$ with \f$ (x,y)\f$ lexicographically.
*/
result_type operator()(
AlgebraicKernel_d_2::Algebraic_real_2 a,
AlgebraicKernel_d_2::Coefficient x,
AlgebraicKernel_d_2::Coefficient y);

/// @}

}; /* end AlgebraicKernel_d_2::CompareXY_2 */

