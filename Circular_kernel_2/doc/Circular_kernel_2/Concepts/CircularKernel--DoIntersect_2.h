
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

Testing whether two curves intersect. 

\refines `Kernel::DoIntersect_2`

\sa `Kernel::do_intersect`

*/

class CircularKernel::DoIntersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
determines if two geometric objects of type `Type1` and `Type2` intersect or not,
for all pairs `Type1` and `Type2`, where the types `Type1` and `Type2` can be any of the
following: 
- `CircularKernel::Line_2`
- `CircularKernel::Line_arc_2`
- `CircularKernel::Circle_2`
- `CircularKernel::Circular_arc_2`
*/ 
bool operator() 
(const Type1 & obj1, const Type2 & obj2); 

/// @}

}; /* end CircularKernel::DoIntersect_2 */

