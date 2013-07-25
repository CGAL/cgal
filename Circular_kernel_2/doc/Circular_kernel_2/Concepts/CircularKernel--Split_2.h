
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::Split_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Splits arc `a` at point `p`, which creates arcs `a1` and `a2`. 
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

