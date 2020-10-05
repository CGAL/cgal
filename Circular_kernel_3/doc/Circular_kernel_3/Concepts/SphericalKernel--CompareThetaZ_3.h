
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept


\sa `SphericalKernel::CompareTheta_3`
\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareXYZ_3`
\sa `SphericalKernel::Equal_3`

*/

class SphericalKernel::CompareThetaZ_3 {
public:

/// \name Operations
/// An object of this type must provide:
/// @{

/*!
Compares `p` and `q` according to the lexicographic ordering on \f$ \theta\f$- and \f$ z\f$-coordinates
in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_z_3_object`.
\pre `p` and `q` lie on the context sphere used by the function `SphericalKernel::compare_theta_z_3_object`, but do not coincide with its poles.

*/
Comparison_result operator()
(const SphericalKernel::Circular_arc_point_3 &p,
const SphericalKernel::Circular_arc_point_3 &q );

/// @}

}; /* end SphericalKernel::CompareThetaZ_3 */

