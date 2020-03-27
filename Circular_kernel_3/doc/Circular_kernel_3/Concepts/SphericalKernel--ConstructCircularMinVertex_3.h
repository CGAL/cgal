
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

*/
class SphericalKernel::ConstructCircularMinVertex_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs the minimal vertex of `l` with lexicographically
smallest coordinates.
*/
SphericalKernel::Circular_arc_point_3 operator()
(const SphericalKernel::Line_arc_3 & l);

/// @}

}; /* end SphericalKernel::ConstructCircularMinVertex_3 */

