
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

An object `fo` of this type must provide: 

*/

class CircularKernel::IsYMonotone_2 {
public:

/// \name Operations 
/// @{

/*! 
Tests whether the arc is \f$ y\f$-monotone. 
*/ 
bool operator() 
(const CircularKernel::Circular_arc_2 & c); 

/*! 
For a line arc, always returns `true`. 
*/ 
bool operator() 
(const CircularKernel::Line_arc_2 & l); 

/// @}

}; /* end CircularKernel::IsYMonotone_2 */

