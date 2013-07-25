
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines `Kernel::HasOnUnboundedSide_2`

\sa `CircularKernel::BoundedSide_2`
\sa `CircularKernel::HasOnBoundedSide_2`

*/

class CircularKernel::HasOnUnboundedSide_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!

*/ 
bool 
operator() 
(const CircularKernel::Circle_2& s, 
const CircularKernel::Circular_arc_point_2& p); 

/// @}

}; /* end CircularKernel::HasOnUnboundedSide_2 */

