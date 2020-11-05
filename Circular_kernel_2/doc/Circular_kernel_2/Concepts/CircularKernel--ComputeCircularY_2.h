
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::ComputeCircularY_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Computes the `y`-coordinate of the point.
*/
CircularKernel::Root_of_2
operator()(const CircularKernel::Circular_arc_point_2 &p);

/// @}

}; /* end CircularKernel::ComputeCircularY_2 */

