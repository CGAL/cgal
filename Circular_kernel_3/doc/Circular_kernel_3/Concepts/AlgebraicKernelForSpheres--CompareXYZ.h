
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

\sa `AlgebraicKernelForSpheres::CompareX`
\sa `AlgebraicKernelForSpheres::CompareY`
\sa `AlgebraicKernelForSpheres::CompareZ`
\sa `AlgebraicKernelForSpheres::CompareXY`

*/

class AlgebraicKernelForSpheres::CompareXYZ {
public:

/// \name Operations
/// A model `fo` of this concept must provide:
/// @{

/*! 
Compares two `Root_for_spheres_2_3` lexicographically. 
*/ 
CGAL::Comparison_result 
operator()(const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r1, 
const AlgebraicKernelForSpheres::Root_for_spheres_2_3 & r2); 

/// @}

}; /* end AlgebraicKernelForSpheres::CompareXYZ */

