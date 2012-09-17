
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalconcept

\refines `Kernel::Intersect_2`

\sa `CGAL::intersection`

*/

class CircularKernel::Intersect_2 {
public:

/// \name Operations
/// A model `fo` of this type must provide: where `Type_1` and
/// `Type_2` can both be either <UL> <LI> `CircularKernel::Line_2` or
/// <LI> `CircularKernel::Line_arc_2` or <LI>
/// `CircularKernel::Circle_2` or <LI>
/// `CircularKernel::Circular_arc_2`. </UL> Depending on the types
/// `Type_1` and `Type_2`, these elements can be assigned to <UL> <LI>
/// `std::pair<CircularKernel::Circular_arc_point_2, unsigned>`, where
/// the unsigned integer is the multiplicity of the corresponding
/// intersection point between `obj_1` and `obj_2`, <LI>
/// `CircularKernel::Circular_arc_2` in case of an overlap of two
/// circular arcs, <LI> `CircularKernel::Line_arc_2` in case of an
/// overlap of two line segments or <LI> `CircularKernel::Line_2` or
/// `CircularKernel::Circle_2` in case of two equal input lines or
/// circles. </UL>
/// @{

/*! 
Copies in the output iterator the intersection elements between the 
two objects. `intersections` iterates on 
elements of type `CGAL::Object`, in lexicographic order. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const Type1 &obj1, const Type2 &obj2, 
OutputIterator intersections); 

/// @}

}; /* end CircularKernel::Intersect_2 */

