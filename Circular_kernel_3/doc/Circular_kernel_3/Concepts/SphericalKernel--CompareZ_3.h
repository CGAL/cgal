
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

\refines `Kernel::CompareZ_3`

\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareXYZ_3`
\sa `SphericalKernel::CompareTheta_3`
\sa `SphericalKernel::CompareThetaZ_3`
\sa `SphericalKernel::Equal_3`

*/

class SphericalKernel::CompareZ_3 {
public:

/// \name Operations
/// An object `fo` of this type must provide in addition:
/// @{

/*! 
Compares the \f$ z\f$-coordinates of \f$ p\f$ and \f$ q\f$. 
*/ 
Comparison_result operator() 
(const SphericalKernel::Circular_arc_point_3 &p, 
const SphericalKernel::Circular_arc_point_3 &q ); 

/// @}

}; /* end SphericalKernel::CompareZ_3 */

