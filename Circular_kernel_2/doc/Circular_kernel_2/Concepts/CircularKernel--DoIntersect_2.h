
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

Testing whether two curves intersect. 

\refines `Kernel::DoIntersect_2`

\sa `Kernel::do_intersect`

*/

class CircularKernel::DoIntersect_2 {
public:

/// \name Operations
/// An object `fo` of this type must provide: for all pairs `Type1`
/// and `Type2`, where the types `Type1` and `Type2` can be any of the
/// following: <UL> <LI> `CircularKernel::Line_2` <LI>
/// `CircularKernel::Line_arc_2` <LI> `CircularKernel::Circle_2` <LI>
/// `CircularKernel::Circular_arc_2` </UL>
/// @{

/*! 
determines if two geometric objects of type Type1 and Type2 intersect or not. 
*/ 
bool operator() 
(const Type1 & obj1, const Type2 & obj2); 

/// @}

}; /* end CircularKernel::DoIntersect_2 */

