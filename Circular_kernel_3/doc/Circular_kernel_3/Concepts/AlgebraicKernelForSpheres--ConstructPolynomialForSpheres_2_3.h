
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

\sa `SphericalKernel::ConstructSphere_3`
\sa `SphericalKernel::GetEquation`

*/

class AlgebraicKernelForSpheres::ConstructPolynomialForSpheres_2_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs polynomial `(x-a)^2 + (y-b)^2 + (z-c)^2 - rsq`.
*/
AlgebraicKernelForSpheres::PolynomialForSpheres_2_3
operator()(const AlgebraicKernelForSpheres::FT &a,
const AlgebraicKernelForSpheres::FT &b,
const AlgebraicKernelForSpheres::FT &c,
const AlgebraicKernelForSpheres::FT &rsq);

/// @}

}; /* end AlgebraicKernelForSpheres::ConstructPolynomialForSpheres_2_3 */

