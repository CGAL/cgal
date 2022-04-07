
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\cgalRefines `Kernel::Intersect_3`
\sa \link intersection_grp `CGAL::intersection()` \endlink

*/
class SphericalKernel::Intersect_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Copies in the output iterator the intersection elements between the
two objects. `intersections` iterates on
elements of type `CGAL::Object`, in lexicographic order
when this ordering is defined on the computed objects.

`Type1` and `Type2` can both be either:

- SphericalKernel::Sphere_3,
- SphericalKernel::Plane_3,
- SphericalKernel::Line_3,
- SphericalKernel::Circle_3,
- SphericalKernel::Line_arc_3 or
- SphericalKernel::Circular_arc_3

depending on the types Type1 and Type2, the computed
`CGAL::Object`'s can be assigned to

- std::pair< \ref SphericalKernel::Circular_arc_point_3, unsigned>,
  where the unsigned integer is the multiplicity of the corresponding
  intersection point between obj_1 and obj_2,
- Type1, when Type1 and Type2 are equal, and
  if the two objects obj1 and obj2 are equal,
- SphericalKernel::Line_3 or SphericalKernel::Circle_3
  when Type1 and Type2 are two-dimensional objects intersecting
  along a curve (2 planes, or 2 spheres, or one plane and one sphere),
- SphericalKernel::Circular_arc_3 in case of an overlap of
  two circular arcs or
- SphericalKernel::Line_arc_3 in case of an overlap of two
  line segments.
*/
template < class OutputIterator >
OutputIterator
operator()(const Type1 &obj1, const Type2 &obj2,
OutputIterator intersections);

/*!
Copies in the output iterator the intersection elements between the
three objects. `intersections` iterates on
elements of type `CGAL::Object`, in lexicographic order
when this ordering is defined on the computed objects.

Type1, Type2 and Type3 can be either:

- SphericalKernel::Sphere_3 or
- SphericalKernel::Plane_3


and depending of these types, the computed `CGAL::Object`'s can be
assigned to

- std::pair< \ref SphericalKernel::Circular_arc_point_3, unsigned>,
  where the unsigned integer is the multiplicity of the corresponding
  intersection point,
- SphericalKernel::Circle_3 or
- Type1, when Type1, Type2 and Type3
  are equal, and if the three objects obj1 and obj2 and obj3
  are equal.
*/
template < class OutputIterator >
OutputIterator
operator()(const Type1 &obj1, const Type2 &obj2,
const Type3 &obj3,
OutputIterator intersections);

///@}


}; /* end SphericalKernel::Intersect_3 */
