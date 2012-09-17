
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class CircularKernel::Split_2 {
public:

/// \name See Also 
/// @{

/*! 
Splits arc \f$ a\f$ at point \f$ p\f$, which creates arcs \f$ a1\f$ and \f$ a2\f$. 
\pre The point lies on the input arc. 
*/ 
void operator() 
(const CircularKernel::Circular_arc_2 &a, 
const CircularKernel::Circular_arc_point_2 &p, 
CircularKernel::Circular_arc_2 &a1, 
CircularKernel::Circular_arc_2 &a2); 

/*! 
Same for a line arc. 
*/ 
void operator() 
(const CircularKernel::Line_arc_2 &l, 
const CircularKernel::Circular_arc_point_2 &p, 
CircularKernel::Line_arc_2 &l1, CircularKernel::Line_arc_2 &l2); 

/// @}

}; /* end CircularKernel::Split_2 */

