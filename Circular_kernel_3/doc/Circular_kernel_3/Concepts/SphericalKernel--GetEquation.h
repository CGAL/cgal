
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept



\sa `SphericalKernel::ConstructPlane_3`
\sa `SphericalKernel::ConstructSphere_3`
\sa `SphericalKernel::ConstructLine_3`
\sa `SphericalKernel::ConstructCircle_3`

*/

class SphericalKernel::GetEquation {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Returns the equation of the plane.
*/
SphericalKernel::Polynomial_1_3
operator()(const SphericalKernel::Plane_3 & p);

/*!
Returns the equation of the sphere.
*/
SphericalKernel::Polynomial_for_spheres_2_3
operator()(const SphericalKernel::Sphere_3 & p);

/*!
Returns the equations of the line.
*/
SphericalKernel::Polynomials_for_line_3
operator()(const SphericalKernel::Line_3 & c);

/*!
Returns the equations of the circle.
*/
SphericalKernel::Polynomials_for_circle_3
operator()(const SphericalKernel::Circle_3 & c);

/// @}

}; /* end SphericalKernel::GetEquation */

