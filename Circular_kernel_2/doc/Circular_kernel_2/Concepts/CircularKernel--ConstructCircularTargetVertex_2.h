
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::ConstructCircularTargetVertex_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{


/*!
Constructs the target vertex of `c`. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Circular_arc_2 & c); 

/*!
Same, for a line segment. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Line_arc_2 & l); 

/// @}

}; /* end CircularKernel::ConstructCircularTargetVertex_2 */

