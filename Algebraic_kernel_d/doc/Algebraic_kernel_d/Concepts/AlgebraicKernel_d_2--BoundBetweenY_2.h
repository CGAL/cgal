
/*!
\ingroup PkgAlgebraicKernelDConceptsBi
\cgalConcept

Computes a number of type
`AlgebraicKernel_d_1::Bound` in-between the second coordinates of two
`AlgebraicKernel_d_2::AlgebraicReal_2`.

\cgalRefines{AdaptableBinaryFunction}

\sa `AlgebraicKernel_d_2::BoundBetweenX_2`

*/

class AlgebraicKernel_d_2::BoundBetweenY_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef AlgebraicKernel_d_1::Bound result_type;

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

Computes a number of type `AlgebraicKernel_d_1::Bound`
in-between the second coordinates of \f$ a\f$ and \f$ b\f$.
\pre \f$ a_y \neq b_y\f$.

*/
result_type
operator()(const first_argument_type & a,
const second_argument_type & b);

/// @}

}; /* end AlgebraicKernel_d_2::BoundBetweenY_2 */

