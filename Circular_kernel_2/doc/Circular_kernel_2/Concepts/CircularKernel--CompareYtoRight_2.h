
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

A function object concept to compare vertically two arcs on the right side of a common point `p`: 

*/

class CircularKernel::CompareYtoRight_2 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{


/*!
For two circular arcs. 
\pre `p` must be a common point to `a1` and `a2`, and `a1` and `a2` must be defined to the right of `p`. 
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

