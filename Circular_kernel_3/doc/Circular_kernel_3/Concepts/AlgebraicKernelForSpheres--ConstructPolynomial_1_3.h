
/*!
\ingroup PkgSphericalKernel3AlgebraicConcepts
\cgalConcept

\sa `SphericalKernel::ConstructPlane_3`
\sa `SphericalKernel::GetEquation`

*/

class AlgebraicKernelForSpheres::ConstructPolynomial_1_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Constructs polynomial `ax+by+cz+d`. 
*/ 
AlgebraicKernelForSpheres::Polynomial_1_3 
operator()(const AlgebraicKernelForSpheres::RT &a, 
const AlgebraicKernelForSpheres::RT &b, 
const AlgebraicKernelForSpheres::RT &c, 
const AlgebraicKernelForSpheres::RT &d); 

/// @}

}; /* end AlgebraicKernelForSpheres::ConstructPolynomial_1_3 */

