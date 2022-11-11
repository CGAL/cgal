
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines `Kernel::ConstructPlane_3`

\sa `SphericalKernel::GetEquation`

*/

class SphericalKernel::ConstructPlane_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs the plane containing the arc.
*/
SphericalKernel::Plane_3 operator()
(const SphericalKernel::Circular_arc_3 &a);

/*!
Constructs a plane from an equation.
*/
SphericalKernel::Plane_3 operator()
(const SphericalKernel::Polynomial_1_3 &p);

/// @}

}; /* end SphericalKernel::ConstructPlane_3 */

