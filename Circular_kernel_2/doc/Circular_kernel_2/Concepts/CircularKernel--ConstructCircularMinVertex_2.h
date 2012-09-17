
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::ConstructCircularMinVertex_2 {
public:

/// \name See Also 
/// @{

/*! 
Constructs the `x`-minimal vertex of `c`. 
\pre The arc `c` is `x`-monotone. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Circular_arc_2 & c); 

/*! 
Same, for a line segment. 
*/ 
CircularKernel::Circular_arc_point_2 operator() 
(const CircularKernel::Line_arc_2 & l); 

/// @}

}; /* end CircularKernel::ConstructCircularMinVertex_2 */

