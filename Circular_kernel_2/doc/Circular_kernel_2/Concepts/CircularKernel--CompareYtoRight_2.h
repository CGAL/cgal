
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

An object `fo` of this type must provide operators that compare vertically 
two arcs on the right side of a common point \f$ p\f$: 

*/

class CircularKernel::CompareYtoRight_2 {
public:

/// \name See Also 
/// @{

/*! 
For two circular arcs. 
\pre \f$ p\f$ must be a common point to \f$ a_1\f$ and \f$ a_2\f$, and \f$ a_1\f$ and \f$ a_2\f$ must be defined to the right of \f$ p\f$. 
*/ 
Comparison_result operator() 
(const Circular_kernel_2::Circular_arc_2 &a1, 
const Circular_kernel_2::Circular_arc_2 &a2, 
const Circular_kernel_2::Circular_arc_point_2 &p); 

/*! 
Same for two segments. 
*/ 
Comparison_result operator() 
(const Circular_kernel_2::Line_arc_2 &a1, 
const Circular_kernel_2::Line_arc_2 &a2, 
const Circular_kernel_2::Circular_arc_point_2 &p); 

/*! 
For a segment and an arc. 
*/ 
Comparison_result operator() 
(const Circular_kernel_2::Line_arc_2 &a1, 
const Circular_kernel_2::Circular_arc_2 &a2, 
const Circular_kernel_2::Circular_arc_point_2 &p); 

/*! 
Same as previous. 
*/ 
Comparison_result operator() 
(const Circular_kernel_2::Circular_arc_2 &a1, 
const Circular_kernel_2::Line_arc_2 &a2, 
const Circular_kernel_2::Circular_arc_point_2 &p); 

/// @}

}; /* end CircularKernel::CompareYtoRight_2 */

