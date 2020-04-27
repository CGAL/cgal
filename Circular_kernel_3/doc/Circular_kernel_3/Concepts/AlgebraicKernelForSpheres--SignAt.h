
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

*/

class AlgebraicKernelForSpheres::SignAt {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Computes the sign of polynomial `p` evaluated at a root `r`.
*/
CGAL::Sign
operator()(const AlgebraicKernelForSpheres::Polynomial_1_3 & p,
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r);

/*!
Same as previous.
*/
CGAL::Sign
operator()(const AlgebraicKernelForSpheres::Polynomial_for_spheres_2_3 & p,
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r);

/// @}

}; /* end AlgebraicKernelForSpheres::SignAt */

