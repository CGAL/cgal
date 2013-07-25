
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines `Kernel::CompareX_2`

\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::CompareXY_2`
\sa `CircularKernel::Equal_2`

*/

class CircularKernel::CompareX_2 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Compares the `x`-coordinates of `p` and `q`. 
*/ 
Comparison_result operator() 
(const CircularKernel::Circular_arc_point_2 &p, 
const CircularKernel::Circular_arc_point_2 &q ); 

/// @}

}; /* end CircularKernel::CompareX_2 */

