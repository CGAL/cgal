
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

\refines `Kernel::CompareX_2`

\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::CompareXY_2`
\sa `CircularKernel::Equal_2`

*/

class CircularKernel::CompareX_2 {
public:

/// \name Operations
/// An object `fo` of this type must provide in addition:
/// @{

/*! 
Compares the \f$ x\f$-coordinates of \f$ p\f$ and \f$ q\f$. 
*/ 
Comparison_result operator() 
(const CircularKernel::Circular_arc_point_2 &p, 
const CircularKernel::Circular_arc_point_2 &q ); 

/// @}

}; /* end CircularKernel::CompareX_2 */

