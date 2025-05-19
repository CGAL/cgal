
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::ComputeApproximateSquaredLength_3}

*/

class SphericalKernel::ComputeApproximateSquaredLength_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Computes an approximation of the squared length of the arc `a`.
*/
double
operator()(const SphericalKernel::Circular_arc_3 & a);

///@}

}; /* end SphericalKernel::ComputeApproximateSquaredLength_3 */

