
/*!
\ingroup PkgSphericalKernel3AlgebraicConcepts
\cgalConcept

\sa `AlgebraicKernelForSpheres::XCriticalPoints`
\sa `AlgebraicKernelForSpheres::ZCriticalPoints`

*/

class AlgebraicKernelForSpheres::YCriticalPoints {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Copies in the output iterator the `y`-critical points of polynomial 
`p`, as objects of type `AlgebraicKernelForSpheres::Root_for_spheres_2_3`. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForSpheres::Polynomial_for_spheres_2_3 &p, 
OutputIterator res); 

/*!
Computes the `i`th `y`-critical point of polynomial `p`. 
*/ 
template < class OutputIterator > 
AlgebraicKernelForSpheres::Root_for_spheres_2_3 
operator()(const AlgebraicKernelForSpheres::Polynomial_for_spheres_2_3 &p, 
bool i); 

/// @}

}; /* end AlgebraicKernelForSpheres::YCriticalPoints */

