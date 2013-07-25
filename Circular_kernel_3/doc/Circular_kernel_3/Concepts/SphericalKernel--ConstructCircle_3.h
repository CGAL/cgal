
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\sa `SphericalKernel::GetEquation`

*/

class SphericalKernel::ConstructCircle_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Constructs the circle containing the arc. 
*/ 
SphericalKernel::Circle_3 operator() 
(const SphericalKernel::Circular_arc_3 &a); 

/*!
Constructs a circle from an equation. 
*/ 
SphericalKernel::Circle_3 operator() 
(const SphericalKernel::Polynomials_for_circles_3 &p); 

/// @}

}; /* end SphericalKernel::ConstructCircle_3 */

