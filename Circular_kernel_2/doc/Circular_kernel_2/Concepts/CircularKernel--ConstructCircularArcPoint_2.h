
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::ConstructCircularArcPoint_2 {
public:

/// \name See Also 
/// @{

/*! 

*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Root_for_circles_2_2 & r); 

/*! 

*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Point_2 & p); 

/// @}

}; /* end CircularKernel::ConstructCircularArcPoint_2 */

