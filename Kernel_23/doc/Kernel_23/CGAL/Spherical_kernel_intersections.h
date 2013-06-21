namespace CGAL {

/*!
\ingroup do_intersect_spherical_grp


\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel for details on a spherical kernel instantiation.

When using a spherical kernel, in addition to the function overloads documented \ref do_intersect_linear_grp "here",
the following function overloads are also available.
*/
/// @{
/*!
checks whether `obj1` and `obj2` intersect.  Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and `obj2`.
Note that while for a polygon we consider the enclosed domain, for an
object of type `Circle_3` or `Sphere_3` only the curve
or the surface is considered. 

`Type1` and `Type2` can be any of
the following:

- `Line_3<SphericalKernel>`
- `Circle_3<SphericalKernel>`
- `Plane_3<SphericalKernel>`
- `Sphere_3<SphericalKernel>`
- `Line_arc_3<SphericalKernel>`
- `Circular_arc_3<SphericalKernel>`

An example illustrating this is presented in
Chapter \ref Chapter_3D_Spherical_Geometry_Kernel.
*/
bool do_intersect(Type1<SphericalKernel> obj1, Type2<SphericalKernel> obj2);

/*!
checks whether `obj1`, `obj2` and `obj3` intersect.

`Type1`, `Type2` and `Type3` can be:

- `Sphere_3<SphericalKernel>`
- `Plane_3<SphericalKernel>`
*/
bool do_intersect(Type1<SphericalKernel> obj1, Type2<SphericalKernel> obj2, Type3<SphericalKernel> obj3);

/// @}

/*!
\ingroup intersection_spherical_grp

\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel for details on a spherical kernel instantiation.

When using a spherical kernel, in addition to the function overloads documented \ref intersection_linear_grp "here",
the following function overloads are also available.

The iterator versions of those functions can be used in conjunction
with `Dispatch_output_iterator`.

Since both the number of intersections, if any, and their types, depend
on the arguments, the function expects an output iterator on
`cpp11::result_of<Kernel::Intersect_3(Type1, Type2)>::%type`,
as presented below.
*/
/// @{

/*!
Copies in the output iterator the intersection elements between the
two objects. `intersections` iterates on
elements of type `result_of< Intersect_3(SphericalType1, SphericalType2) >`, in lexicographic order,
when this ordering is defined on the computed objects,

where `SphericalType1` and `SphericalType2` can both be one of:

- `Sphere_3<SphericalKernel>`,
- `Plane_3<SphericalKernel>`,
- `Line_3<SphericalKernel>`,
- `Circle_3<SphericalKernel>`,
- `Line_arc_3<SphericalKernel>` or
- `Circular_arc_3<SphericalKernel>`,


and depending on the types `SphericalType1` and `SphericalType2`, the computed 
type can be

- `std::pair<Circular_arc_point_3<SphericalKernel>, unsigned>`,
  where the unsigned integer is the multiplicity of the corresponding
  intersection point between `obj1` and `obj2`,
- `SphericalType1`, when `SphericalType1` and `SphericalType2` are equal, 
  and if the two objets `obj1` and `obj2` are equal,
- `Line_3<SphericalKernel>` or 
  `Circle_3<SphericalKernel>` when `SphericalType1` and `SphericalType2` 
  are two-dimensional objets intersecting along a curve (2 planes, or 2
  spheres, or one plane and one sphere),
- `Circular_arc_3<SphericalKernel>` in case of an overlap of 
  two circular arcs or
- `Line_arc_3<SphericalKernel>` in case of an overlap of two 
  line segments. 
*/
template < typename SphericalType1, typename SphericalType1,  typename OutputIterator >
OutputIterator
intersection(const SphericalType1 &obj1, const SphericalType2 &obj2,
             OutputIterator intersections);


/*!
Copies in the output iterator the intersection elements between the
three objects. `intersections` iterates on
elements of type `boost::variant< Circle_3, Plane_3, Sphere_3, std::pair< Circular_arc_point_3, unsigned > >`, in lexicographic order 
when this ordering is defined on the computed objects

where `Type1`, `Type2` and `Type3`
can be either

- `Sphere_3<SphericalKernel>` or
- `Plane_3<SphericalKernel>`


and depending of these types, the computed return value

- `std::pair<Circular_arc_point_3<SphericalKernel>, unsigned>`,
  where the unsigned integer is the multiplicity of the corresponding
  intersection point,
- `Circle_3<SphericalKernel>` or
- `Type1`, when `Type1`, `Type2` and 
  `Type3` are equal, and if the three objets `obj1` and `obj2` 
  and `obj3` are equal.
*/
  template < typename Type1, typename Type2, typename Type3, typename OutputIterator >
OutputIterator
intersection(const Type1 &obj1, const Type2 &obj2, const Type3 &obj3,
             OutputIterator intersections);

/// @}

} // CGAL


