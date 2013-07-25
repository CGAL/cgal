
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines `Kernel::CompareXY_3`

\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXYZ_3`
\sa `SphericalKernel::CompareTheta_3`
\sa `SphericalKernel::CompareThetaZ_3`
\sa `SphericalKernel::Equal_3`

*/

class SphericalKernel::CompareXY_3 {
public:

/// \name Operations 
///  An object of this type must provide in addition:
/// @{

/*!
Compares `p` and `q` according to the lexicographic ordering on \f$ x\f$- and \f$ y\f$-coordinates. 
*/ 
Comparison_result operator() 
(const SphericalKernel::Circular_arc_point_3 &p, 
const SphericalKernel::Circular_arc_point_3 &q ); 

/// @}

}; /* end SphericalKernel::CompareXY_3 */

