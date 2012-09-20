
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalconcept

\brief

\refines `Kernel::ConstructSphere_3`

\sa `SphericalKernel::GetEquation`

*/

class SphericalKernel::ConstructSphere_3 {
public:

/// \name Operations 
/// A model `fo` of this concept must provide:
/// @{

/*! 
Returns the diametral sphere of the supporting circle of arc `a`. 
*/ 
SphericalKernel::Sphere_3 
operator()(const SphericalKernel::Circular_arc_3 & a); 

/*! 
Constructs a sphere from an equation. 
*/ 
SphericalKernel::Sphere_3 operator() 
(const SphericalKernel::Polynomial_2_3 &p); 

/// @}

}; /* end SphericalKernel::ConstructSphere_3 */

