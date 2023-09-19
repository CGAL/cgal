
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::CompareXYZ_3}

\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareTheta_3`
\sa `SphericalKernel::CompareThetaZ_3`
\sa `SphericalKernel::Equal_3`

*/
class SphericalKernel::CompareXYZ_3 {
public:

/// \name Operations
/// An object of this type must provide in addition:
/// @{

/*!
Compares `p` and `q` according to the lexicographic ordering on \f$ x\f$-, \f$ y\f$-,
and \f$ z\f$-coordinates.
*/
Comparison_result operator()
(const SphericalKernel::Circular_arc_point_3 &p,
const SphericalKernel::Circular_arc_point_3 &q );

/// @}

}; /* end SphericalKernel::CompareXYZ_3 */

