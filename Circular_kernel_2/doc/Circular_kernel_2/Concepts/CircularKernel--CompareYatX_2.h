
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

An object `fo` of this type must provide two operators that compare 
a point \f$ p\f$ and an arc \f$ a\f$ on the vertical line passing through \f$ p\f$. 

*/

class CircularKernel::CompareYatX_2 {
public:

/// \name See Also 
/// @{

/*! 
For a circular arc. 
\pre The arc \f$ a\f$ must be monotone and \f$ p\f$ must be in the vertical range of \f$ a\f$. 
*/ 
Comparison_result operator() 
(const CircularKernel::Circular_arc_point_2 &p, 
const CircularKernel::Circular_arc_2 &a); 

/*! 
Same for a segment. 
*/ 
Comparison_result operator() 
(const CircularKernel::Circular_arc_point_2 &p, 
const CircularKernel::Line_arc_2 &a); 

/// @}

}; /* end CircularKernel::CompareYatX_2 */

