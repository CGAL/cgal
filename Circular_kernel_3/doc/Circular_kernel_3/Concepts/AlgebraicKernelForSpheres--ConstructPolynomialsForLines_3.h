
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

\sa `SphericalKernel::ConstructLine_3`
\sa `SphericalKernel::GetEquation`

*/

class AlgebraicKernelForSpheres::ConstructPolynomialsForLines_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs polynomials `a1.t+a2`, `b1.t+b2`, `c1.t+c2`.
*/
AlgebraicKernelForSpheres::Polynomials_for_line_3
operator()(const AlgebraicKernelForSpheres::FT &a1,
const AlgebraicKernelForSpheres::FT &a2,
const AlgebraicKernelForSpheres::FT &b1,
const AlgebraicKernelForSpheres::FT &b2,
const AlgebraicKernelForSpheres::FT &c1,
const AlgebraicKernelForSpheres::FT &c2);

/// @}

}; /* end AlgebraicKernelForSpheres::ConstructPolynomialsForLines_3 */

