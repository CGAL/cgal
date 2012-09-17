
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

Testing equality between objects. 

\refines `Kernel::Equal_2`

\sa `CircularKernel::CompareX_2`
\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::CompareXY_2`

*/

class CircularKernel::Equal_2 {
public:

/// \name Operations
/// An object `fo` of this type must provide in addition: For the
/// sake of completeness, the `operator()` must also be defined for a
/// `Line_arc_2` and a `Circular_arc_2` as arguments (in any order),
/// and it always returns `false`.
/// @{

/*! 
For two points. 
*/ 
bool operator() 
(const CircularKernel::Circular_arc_point_2 &p0, 
const CircularKernel::Circular_arc_point_2 &p1); 

/*! 
For two arcs. 
*/ 
bool operator() 
(const CircularKernel::Circular_arc_2 &a0, 
const CircularKernel::Circular_arc_2 &a1); 

/*! 
For two segments. 
*/ 
bool operator() 
(const CircularKernel::Line_arc_2 &a0, 
const CircularKernel::Line_arc_2 &a1); 

/// @}

}; /* end CircularKernel::Equal_2 */

