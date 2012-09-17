
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::ConstructCircularMaxVertex_2 {
public:

/// \name See Also 
/// @{

/*! 
Constructs the \f$ x\f$-maximal vertex of `c`. 
\pre The arc `c` is \f$ x\f$-monotone. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Circular_arc_2 & c); 

/*! 
Same, for a line segment. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Line_arc_2 & l); 

/// @}

}; /* end CircularKernel::ConstructCircularMaxVertex_2 */

