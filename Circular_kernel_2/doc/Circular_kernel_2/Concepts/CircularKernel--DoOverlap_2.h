
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

Testing whether the interiors of two curves overlap. 

\refines ::Kernel::DoOverlap_2 

*/

class CircularKernel::DoOverlap_2 {
public:

/// \name Operations 
/// An object `fo` of this type must provide: 
/// @{

/*! 
For two line arcs. 
*/ 
bool operator() 
(const CircularKernel::Line_arc_2 & l0, 
const CircularKernel::Line_arc_2 & l1); 

/*! 
For two circular arcs. 
*/ 
bool operator() 
(const CircularKernel::Circular_arc_2 & a0, 
const CircularKernel::Circular_arc_2 & a1); 

/// @}

}; /* end CircularKernel::DoOverlap_2 */

