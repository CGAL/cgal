
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

\refines `Kernel::CompareXY_2`

\sa `CircularKernel::CompareX_2`
\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::Equal_2`

*/

class CircularKernel::CompareXY_2 {
public:

/// \name Operations 
/// An object `fo` of this type must provide in addition:
/// @{

/*! 
Compares \f$ p\f$ and \f$ q\f$ according to the lexicographic ordering on \f$ x\f$- and \f$ y\f$-coordinates. 
*/ 
Comparison_result operator() 
(const CircularKernel::Circular_arc_point_2 &p, 
const CircularKernel::Circular_arc_point_2 &q ); 

/// @}

}; /* end CircularKernel::CompareXY_2 */

