
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForSpheres::CompareX`
\sa `AlgebraicKernelForSpheres::CompareY`
\sa `AlgebraicKernelForSpheres::CompareXY`
\sa `AlgebraicKernelForSpheres::CompareXYZ`

*/

class AlgebraicKernelForSpheres::CompareZ {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the \f$ z\f$-coordinates of two `Root_for_spheres_2_3`.
*/
CGAL::Comparison_result
operator()(const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r1,
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r2);

/// @}

}; /* end AlgebraicKernelForSpheres::CompareZ */

