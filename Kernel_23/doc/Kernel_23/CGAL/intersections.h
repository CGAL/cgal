namespace CGAL {

/*!
\addtogroup do_intersect

\brief 
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.
*/

/*!
\addtogroup do_intersect_linear do_intersect (2D/3D Linear Kernel)
\ingroup do_intersect

\code
#include <CGAL/intersections.h>
\endcode
  
\sa \ref do_intersect_circular
\sa \ref do_intersect_spherical
\sa `intersection`
  
\details See Chapter  \ref chapterkernel23 for details on a linear kernel instantiation.
*/
/// @{
/*!
checks whether `obj1` and `obj2` intersect.  Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and `obj2`.
Note that for objects like triangles and polygons that enclose a
bounded region, this region is part of the object.

The types `Type1` and `Type2` can be any of the following:

- `Point_2<Kernel>`
- `Line_2<Kernel>`
- `Ray_2<Kernel>`
- `Segment_2<Kernel>`
- `Triangle_2<Kernel>`
- `Iso_rectangle_2<Kernel>`

Also, `Type1` and `Type2` can be both of type

- `Line_2<Kernel>`
- `Circle_2<Kernel>`

In three-dimensional space, the types `Type1` and
`Type2` can be any of the following:

- `Plane_3<Kernel>`
- `Line_3<Kernel>`
- `Ray_3<Kernel>`
- `Segment_3<Kernel>`
- `Triangle_3<Kernel>`.
- `Bbox_3`.

Also, `Type1` and `Type2` can be respectively of types

- `Triangle_3<Kernel>` and `Tetrahedron_3<Kernel>`
- `Plane_3<Kernel>` and `Sphere_3<Kernel>` (or the contrary)
- `Sphere_3<Kernel>` and `Sphere_3<Kernel>`.

*/
bool do_intersect(Type1<Kernel> obj1, Type2<Kernel> obj2);
/// @}



/*!
\addtogroup do_intersect_circular do_intersect (2D Circular Kernel)
\ingroup do_intersect

\code
#include <CGAL/Circular_kernel_intersections.h>
\endcode
  
\sa \ref do_intersect_linear
\sa \ref do_intersect_spherical
\sa `intersection`

\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel for details on a circular kernel instantiation.


When using a circular kernel, in addition to the function overloads documented \ref do_intersect_linear "here",
the following function overloads are also available.


*/
/// @{
/*!
checks whether `obj1` and `obj2` intersect.  Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and `obj2`.
Note that for objects like triangles and polygons that enclose a
bounded region, this region is part of the object.

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
\addtogroup do_intersect_spherical do_intersect (3D Spherical Kernel)
\ingroup do_intersect

\code
#include <CGAL/Spherical_kernel_intersections.h>
\endcode
  
\sa \ref do_intersect_linear
\sa \ref do_intersect_circular
\sa `intersection`

\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel for details on a spherical kernel instantiation.


When using a circular kernel, in addition to the function overloads documented \ref do_intersect_linear "here",
the following function overloads are also available.


*/
/// @{
/*!
checks whether `obj1` and `obj2` intersect.  Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and `obj2`.
Note that for objects like triangles and polygons that enclose a
bounded region, this region is part of the object.

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
\addtogroup intersection

The macro `CGAL_INTERSECTION_VERSION` can be used to configure
which version of the `intersection()` function should be used and
enables the corresponding APIs in other \cgal packages. It should be
defined before any \cgal header is included.

- `CGAL_INTERSECTION_VERSION == 1` `intersection()` uses `Object`
- `CGAL_INTERSECTION_VERSION == 2` `intersection()` uses `boost::optional< boost::variant< T... > >`

*/

#define CGAL_INTERSECTION_VERSION


/*!
\addtogroup intersection

\brief 
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.

### Notes on Backward Compatibility ###

The `intersection()` function used to return an `Object`, but starting with 
\cgal 4.2 the
return type is determined by a metafunction defined by the kernel. To
preserve backward compatibility `Object` can be
constructed from the new return types implicitly, but switching to the
new style is recommended. To enable the old style without any overhead,
the macro `CGAL_INTERSECTION_VERSION` must be defined to
`1` before any \cgal header is included.

### Upgrading code from CGAL::Object to boost::variant###

Code can be upgraded by using either `boost::get` or the
`boost::static_visitor<T>`.

\code
#include <CGAL/intersections.h>

template<typename R>
struct Intersection_visitor : public boost::static_visitor<> {
  void operator()(const Point_2& p) const 
  { // Point_2  
  }
  void operator()(const Segment_2& s) const 
  { // Segment_2 
  }
};

template <class R>
void foo(Segment_2<R> seg, Line_2<R> lin)
{
  CGAL::Object obj = intersection(seg1, seg2);
  if(const Point_2* foo = object_cast<Point_2>(&obj)) {
    // Point_2 
  } else if(const Segment_2* foo = object_cast<Segment_2>(&obj)) {
    // Segment_2 
  } else {
    // empty
  }

  // becomes
  auto result = intersection(seg, lin);
  if(result) {
    if(const Point_2* foo = boost::get<Point_2>(&*obj)) {
      // Point_2 
    } else if(const Segment_2* foo = boost::get<Segment_2>(&*obj)) {
      // Segment_2 
    }
  } else {
    // empty
  }

  // or with boost::static_visitor<T>
  if(result) { boost::apply_visitor(Intersection_visitor(), *result); } 
  else { // empty  
  }
}
\endcode



*/


/*!
\addtogroup intersection_linear intersection (2D/3D Linear Kernel)
\ingroup intersection

\code
#include <CGAL/intersections.h>
\endcode

\sa intersection_circular
\sa intersection_spherical
\sa `CGAL::do_intersect` 

\details See Chapter  \ref chapterkernel23 for details on a linear kernel instantiation.

*/
/// @{

/*!
Two objects `obj1` and `obj2` intersect if there is a point `p` that
is part of both `obj1` and `obj2`.  The intersection region of those
two objects is defined as the set of all points `p` that are part of
both `obj1` and `obj2`.  Note that for objects like triangles and
polygons that enclose a bounded region, this region is considered part
of the object.  If a segment lies completely inside a triangle, then
those two objects intersect and the intersection region is the
complete segment. 

Here, `Intersect_23` means either `Intersect_2` or `Intersect_3`,
depending on the arguments.

The possible values for types `Type1` and `Type2` and
the value for `T...` in `boost::optional< boost::variant<
  T... > >` are the following and can be obtained through
`boost::result_of(Intersect_2(A, B)>::type` or
`boost::result_of(Intersect_3(A, B)>::type`.
 


<DIV ALIGN="CENTER">
<TABLE CELLPADDING=3 BORDER="1">
<TR> <TH> Type1 </TH>
 <TH> Type2 </TH>
 <TH> Return Type:  `T...` </TH>
</TR>
<TR>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD>Iso_rectangle_2</TD>
</TR>

<TR>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Iso_rectangle_2 </TD>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD>Point_2, or Segment_2, or Triangle_2, or std::vector&lt;Point_2&gt;</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD>Point_2, or Line_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD>Point_2, or Ray_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_2 </TD>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD>Point_2, or Segment_2, or Ray_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_2 </TD>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Segment_2 </TD>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD>Point_2, or Segment_2</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD VALIGN="CENTER" > Triangle_2 </TD>
    <TD>Point_2, or Segment_2, or Triangle_2, or std::vector&lt;Point_2&gt;</TD>
</TR>

</TABLE>
<TABLE CELLPADDING=3 BORDER="1">
<TR> <TH> Type1 </TH>
 <TH> Type2 </TH>
 <TH> Return Type: `T...` </TH>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD>Point_3, or Line_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD>Point_3, or Line_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD>Point_3, or Ray_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD>Line_3, or Plane_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD>Point_3, or Ray_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Sphere_3 </TD>
    <TD>Point_3, or Circle_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD>Point_3, or Ray_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
p    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Sphere_3 </TD>
    <TD VALIGN="CENTER" > Sphere_3 </TD>
    <TD>Point_3, or Circle_3, or Sphere_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3, or std::vector &lt; Point_3  &gt;</TD>
</TR>
</TABLE>
</DIV>

### Example ###

The following example demonstrates the most common use of 
`intersection` routines with the 2D and 3D Linear %Kernel.

\code
#include <CGAL/intersections.h>

template<typename R>
struct Intersection_visitor {
  typedef result_type void;
  void operator()(const Point_2<R>& p) const 
  { // handle point
  }
  void operator()(const Segment_2<R>& s) const 
  { 
    // handle segment 
  }
};

template <class R>
void foo(const Segment_2<R>& seg, const Line_2<R>& lin)
{
  // with C++11 support
  // auto result = intersection(seg, lin);

  // without C++11
  boost::result_of<R::Intersect_2(Segment_2<R>, Line_2<R>)>::type 
    result = intersection(seg, lin);

  if(result) { boost::apply_visitor(Intersection_visitor(), *result); } 
  else { 
    // no intersection 
  }

  // alternatively:
  if(result) {
    if(const Segment_2<R>* s = boost::get<Segment_2>(&*result)) {
      // handle segment
    } else {
      const Point_2<R>* p = boost::get<Point_2>(&*result);
      // handle point
    }
}
\endcode


Another example showing the use of the intersection function as a plain function call and with `Dispatch_output_iterator` combined with a standard library algorithm.

\cgalExample{Kernel_23/intersections.cpp}



*/
boost::result_of<Intersect_23(Type1, Type2)>::type
intersection(Type1<Kernel> obj1, Type2<Kernel> obj2);

/*!
returns the intersection of 3 planes, which can be either a
point, a line, a plane, or empty.
*/
boost::optional< boost::variant< Point_3, Line_3, Plane_3 > > 
intersection(const Plane_3<Kernel>& pl1,
             const Plane_3<Kernel>& pl2,
             const Plane_3<Kernel>& pl3);

/// @}

/*!
\addtogroup intersection_circular intersection (2D Circular Kernel)
\ingroup intersection

\code
#include <CGAL/Circular_kernel_intersections.h>
\endcode

\sa intersection_linear
\sa intersection_spherical
\sa `CGAL::do_intersect` 

\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel for details on a circular kernel instantiation.

When using a circular kernel, in addition to the function overloads documented \ref intersection_linear "here",
the following function overloads are also available.

The iterator versions of those functions can be used in conjunction
with `Dispatch_output_iterator`.

Since both the number of intersections, if any, and their type, depend
on the arguments, the function expects an output iterator on
`boost::result_of<K::Intersect_2(Type1, Type2)>::type`, as
presented below.

*/

/// @{

/*!
Constructs the intersection elements between the
two objects and stores them in the OutputIterator in lexicographic order,

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

/*!
\addtogroup intersection_spherical intersection (3D Spherical Kernel)
\ingroup intersection

\code
#include <CGAL/Spherical_kernel_intersections.h>
\endcode

\sa intersection_linear
\sa intersection_circular
\sa `CGAL::do_intersect` 

\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel for details on a spherical kernel instantiation.

When using a spherical kernel, in addition to the function overloads documented \ref intersection_linear "here",
the following function overloads are also available.

The iterator versions of those functions can be used in conjunction
`Dispatch_output_iterator`.

Since both the number of intersections, if any, and their type, depend
on the arguments, the functions expects an output iterator on
`boost::result_of<Kernel::Intersect_3(Type1, Type2)>::type`,
as presented below.
*/
/// @{

/*!
Copies in the output iterator the intersection elements between the
two objects. `intersections` iterates on
elements of type `result_of< Intersect_3(SphericalType1, SphericalType2) >`, in lexicographic order,
when this ordering is defined on the computed objects,

where `SphericalType1` and `SphericalType2` can both be either:

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
}
