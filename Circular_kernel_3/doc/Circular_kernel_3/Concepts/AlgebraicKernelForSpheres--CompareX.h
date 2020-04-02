
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForSpheres::CompareY`
\sa `AlgebraicKernelForSpheres::CompareZ`
\sa `AlgebraicKernelForSpheres::CompareXY`
\sa `AlgebraicKernelForSpheres::CompareXYZ`

*/

class AlgebraicKernelForSpheres::CompareX {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the `x` (first) variables of two `Root_for_spheres_2_3`.
*/
CGAL::Comparison_result
operator()(const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r1,
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r2);

/// @}

}; /* end AlgebraicKernelForSpheres::CompareX */

