
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines{Kernel::CompareY_3}

\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareXYZ_3`
\sa `SphericalKernel::CompareTheta_3`
\sa `SphericalKernel::CompareThetaZ_3`
\sa `SphericalKernel::Equal_3`

*/

class SphericalKernel::CompareY_3 {
public:

/// \name Operations
/// An object of this type must provide in addition:
/// @{

/*!
Compares the \f$ y\f$-coordinates of `p` and `q`.
*/
Comparison_result operator()
(const SphericalKernel::Circular_arc_point_3 &p,
const SphericalKernel::Circular_arc_point_3 &q );

/// @}

}; /* end SphericalKernel::CompareY_3 */

