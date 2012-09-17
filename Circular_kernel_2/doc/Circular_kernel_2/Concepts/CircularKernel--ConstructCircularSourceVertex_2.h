
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::ConstructCircularSourceVertex_2 {
public:

/// \name See Also 
/// @{

/*! 
Constructs the source vertex of `c`. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Circular_arc_2 & c); 

/*! 
Same, for a line segment. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Line_arc_2 & l); 

/// @}

}; /* end CircularKernel::ConstructCircularSourceVertex_2 */

