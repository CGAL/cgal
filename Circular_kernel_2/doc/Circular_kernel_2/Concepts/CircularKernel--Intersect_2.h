
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines `Kernel::Intersect_2`

\sa \link intersection_grp `CGAL::intersection()` \endlink

*/

class CircularKernel::Intersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Copies in the output iterator the intersection elements between the
two objects. `intersections` iterates on
elements of type `CGAL::Object`, in lexicographic order.

`Type_1` and `Type_2` can both be either:
- `CircularKernel::Line_2`
- `CircularKernel::Line_arc_2`
- `CircularKernel::Circle_2`
- `CircularKernel::Circular_arc_2`.

Depending on the types `Type_1` and `Type_2`, these elements can be assigned to
- `std::pair<CircularKernel::Circular_arc_point_2, unsigned>`, where the unsigned integer is the multiplicity of the corresponding intersection point between `obj_1` and `obj_2`,
- `CircularKernel::Circular_arc_2` in case of an overlap of two circular arcs,
- `CircularKernel::Line_arc_2` in case of an overlap of two line segments or
- `CircularKernel::Line_2` or `CircularKernel::Circle_2` in case of two equal input lines or circles.
*/
template < class OutputIterator >
OutputIterator
operator()(const Type1 &obj1, const Type2 &obj2,
OutputIterator intersections);

/// @}

}; /* end CircularKernel::Intersect_2 */

