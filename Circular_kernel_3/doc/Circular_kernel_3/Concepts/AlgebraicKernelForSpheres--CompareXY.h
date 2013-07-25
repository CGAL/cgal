
/*!
\ingroup PkgSphericalKernel3AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForSpheres::CompareX`
\sa `AlgebraicKernelForSpheres::CompareY`
\sa `AlgebraicKernelForSpheres::CompareZ`
\sa `AlgebraicKernelForSpheres::CompareXYZ`

*/

class AlgebraicKernelForSpheres::CompareXY {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the `x` and `y` (the two first) variables of 
two `Root_for_spheres_2_3` lexicographically. 
*/ 
CGAL::Comparison_result 
operator()(const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r1, 
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r2); 

/// @}

}; /* end AlgebraicKernelForSpheres::CompareXY */

