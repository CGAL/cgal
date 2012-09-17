
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::ComputeCircularX_2 {
public:

/// \name See Also 
/// @{

/*! 
Computes the \f$ x\f$-coordinate of the point. 
*/ 
CircularKernel::Root_of_2 
operator()(const CircularKernel::Circular_arc_point_2 &p); 

/// @}

}; /* end CircularKernel::ComputeCircularX_2 */

