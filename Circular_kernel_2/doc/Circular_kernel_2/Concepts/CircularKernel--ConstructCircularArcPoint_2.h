
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::ConstructCircularArcPoint_2 {
public:

/// \name Operations
/// A model of this concept must provide:
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

