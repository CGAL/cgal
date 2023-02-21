namespace CGAL {

/*!
\ingroup do_intersect_circular_grp


\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel for details on a circular kernel instantiation.

When using a circular kernel, in addition to the function overloads documented \ref do_intersect_linear_grp "here",
the following function overloads are also available.
*/

/// @{

/*!
checks whether `obj1` and `obj2` intersect.  Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and `obj2`.
Note that while for a polygon we consider the enclosed domain, for an
object of type `Circle_2` only the curve is considered.

`Type1` and `Type2` can be any of
the following:

- `Line_2<CircularKernel>`
- `Circle_2<CircularKernel>`
- `Line_arc_2<CircularKernel>`
- `Circular_arc_2<CircularKernel>`

An example illustrating this is presented in
Chapter  \ref Chapter_2D_Circular_Geometry_Kernel.
*/
bool do_intersect(Type1<CircularKernel> obj1, Type2<CircularKernel> obj2);
/// @}

/*!
\ingroup intersection_circular_grp


\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel for details on a circular kernel instantiation.

When using a circular kernel, in addition to the function overloads documented \ref intersection_linear_grp "here",
the following function overloads are also available.

The iterator versions of those functions can be used in conjunction
with `Dispatch_output_iterator`.

Since both the number of intersections, if any, and types of the intersection results
depend on the arguments, the function expects an output iterator on `K::Intersect_2(Type1, Type2)`
as presented below.
*/

/// @{

/*!

Constructs the intersection elements between the two input
objects and stores them in the OutputIterator in lexicographic order,
where both, `Type1` and `Type2`, can be either

- `Line_2<CircularKernel>` or
- `Line_arc_2<CircularKernel>` or
- `Circle_2<CircularKernel>` or
- `Circular_arc_2<CircularKernel>`

\details Depending on the types `Type1` and `Type2`, these elements can be assigned to

- `std::pair<Circular_arc_point_2<CircularKernel>, unsigned>`,
  where the unsigned integer is the multiplicity of the corresponding
  intersection point between `obj1` and `obj2`,
- `Circular_arc_2<CircularKernel>` in case of an overlap of
  two circular arcs,
- `Line_arc_2<CircularKernel>` in case of an overlap of two
  line segments or
- `Line_2<CircularKernel>` or
  `Circle_2<CircularKernel>` in case of two equal input lines or circles.

*/
template < typename Type1, typename Type2, typename OutputIterator >
OutputIterator
intersection(const Type1 &obj1, const Type2 &obj2,
             OutputIterator intersections);

/// @}

} // CGAL

