/// \file intersections.h

namespace CGAL {

/*!
\addtogroup do_intersect_grp

\brief
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.
*/

/*!
\addtogroup do_intersect_linear_grp
\ingroup do_intersect_grp

\sa `do_intersect_circular_grp`
\sa `do_intersect_spherical_grp`
\sa `intersection_grp`

\details See Chapter  \ref chapterkernel23 "2D and 3D Geometry Kernel" for details on a linear kernel instantiation.
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

- `Bbox_3`.
- `Point_3<Kernel>`
- `Plane_3<Kernel>`
- `Line_3<Kernel>`
- `Ray_3<Kernel>`
- `Segment_3<Kernel>`
- `Sphere_3<Kernel>`
- `Triangle_3<Kernel>`.
- `Tetrahedron_3<Kernel>`.
*/
bool do_intersect(Type1<Kernel> obj1, Type2<Kernel> obj2);

/*!
checks whether `obj1`, `obj2` and `obj3` intersect.
*/
bool do_intersect(Plane_3<Kernel> obj1, Plane_3<Kernel> obj2, Plane_3<Kernel> obj3);

/// @}



/*!
\addtogroup intersection_grp

\brief
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.

*/

/*!
\addtogroup intersection_linear_grp
\ingroup intersection_grp

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

The following tables give the possible values for `Type1` and `Type2`.

\cgalHeading{2D Intersections}

The return type of intersecting two objects of the types `Type1` and `Type2` can be
specified through the placeholder type specifier `auto`. It is equivalent to
`std::optional< std::variant< T... > >`, the last column in the table providing
the template parameter pack.

<DIV ALIGN="CENTER">
<TABLE CELLPADDING=3 BORDER="1">
<TR>
 <TH> Type1 </TH>
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
</DIV>

Additional overloads are provided for the type `Point_2` combined with any other type with the result type being
`std::optional< std::variant< Point_2 > >`.
Overloads are also provided for the type `Bbox_2`, for all
intersections existing with the type `Iso_rectangle_2`. Note that the return type for `Bbox_2` - `Bbox_2`
 is `Bbox_2` and not `Iso_rectangle_2`.

\cgalHeading{3D Intersections}

The return type of intersecting two objects of the types `Type1` and `Type2` can be
specified through the placeholder type specifier `auto`. It is equivalent to
`std::optional< std::variant< T... > >`, the last column in the table providing
the template parameter pack.

<DIV ALIGN="CENTER">
<TABLE CELLPADDING=3 BORDER="1">
<TR>
 <TH> Type1 </TH>
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
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Tetrahedron_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Line_3 </TD>
    <TD VALIGN="CENTER" > Iso_cuboid_3 </TD>
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
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Tetrahedron_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3, or std::vector<Point_3></TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Plane_3 </TD>
    <TD VALIGN="CENTER" > Iso_cuboid_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3, or std::vector<Point_3></TD>
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
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD VALIGN="CENTER" > Tetrahedron_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Ray_3 </TD>
    <TD VALIGN="CENTER" > Iso_cuboid_3 </TD>
    <TD>Point_3, or Segment_3</TD>
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
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD VALIGN="CENTER" > Tetrahedron_3 </TD>
    <TD>Point_3, or Segment_3</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Segment_3 </TD>
    <TD VALIGN="CENTER" > Iso_cuboid_3 </TD>
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
<TR>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD VALIGN="CENTER" > Tetrahedron_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3, or std::vector &lt; Point_3  &gt;</TD>
</TR>
<TR>
    <TD VALIGN="CENTER" > Triangle_3 </TD>
    <TD VALIGN="CENTER" > Iso_cuboid_3 </TD>
    <TD>Point_3, or Segment_3, or Triangle_3, or std::vector &lt; Point_3  &gt;</TD>
</TR>
</TABLE>
</DIV>

Additional overloads are provided for the type `Point_3` combined with any other type with the result type being
`std::optional< std::variant< Point_3 > >`. Overloads are also provided for the type `Bbox_3`, for all
intersections existing with the type `Iso_cuboid_3`. Note that the return type for `Bbox_3` - `Bbox_3`
 is `Bbox_3` and not `Iso_cuboid_3`.


\cgalHeading{Examples}

The following examples demonstrate the most common use of
`intersection()` functions with the 2D and 3D Linear %Kernel.

In the first two examples we intersect a segment and a line.
The result type can be specified through the placeholder type specifier `auto`,
but you must anyway know that the result type is a `std::optional<std::variant<..> >`,
in order to unpack the point or segment.

<A HREF="https://www.boost.org/libs/optional/">`std::optional`</A> comes in
as there might be no intersection. <A HREF="https://www.boost.org/libs/variant/">`std::variant`</A> comes in
as, if there is an intersection, it is either a point or a segment.

As explained in the boost manual pages for <A HREF="https://www.boost.org/libs/variant/">`std::variant`</A>, there are two ways to access the variants. The first examples uses `boost::get`.

\cgalExample{Kernel_23/intersection_get.cpp}

The second example uses `std::visit`.

\cgalExample{Kernel_23/intersection_visitor.cpp}

A third example shows the use of the intersection function as a
plain function call and with `Dispatch_output_iterator`, combined with
a standard library algorithm.

\cgalExample{Kernel_23/intersections.cpp}

*/
template <typename Kernel>
decltype(auto)
intersection(Type1<Kernel> obj1, Type2<Kernel> obj2);

/*!
returns the intersection of 3 planes, which can be a
point, a line, a plane, or empty.
*/
template <typename Kernel>
decltype(auto)
intersection(const Plane_3<Kernel>& pl1,
             const Plane_3<Kernel>& pl2,
             const Plane_3<Kernel>& pl3);

/// @}
}
