/// \file intersections.h

/*!
\def CGAL_INTERSECTION_VERSION

\ingroup intersection_grp

The macro `CGAL_INTERSECTION_VERSION` can be used to configure
which version of the \ref intersection_grp function should be used and
enables the corresponding APIs in other \cgal packages. It should be
defined before any \cgal header is included.

- `CGAL_INTERSECTION_VERSION == 1` \ref intersection_grp uses `Object`
- `CGAL_INTERSECTION_VERSION == 2` \ref intersection_grp uses `boost::optional< boost::variant< T... > >`

*/
#define CGAL_INTERSECTION_VERSION

namespace CGAL {

/*!
\addtogroup do_intersect_grp

\brief
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.
*/

/*!
\addtogroup do_intersect_linear_grp
\ingroup do_intersect

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

- `Plane_3<Kernel>`
- `Line_3<Kernel>`
- `Ray_3<Kernel>`
- `Segment_3<Kernel>`
- `Triangle_3<Kernel>`.
- `Bbox_3`.

Also, `Type1` and `Type2` can be respectively of types

- `Triangle_3<Kernel>` and `Tetrahedron_3<Kernel>`
- `Plane_3<Kernel>` and `Sphere_3<Kernel>` (or the contrary)
- `Sphere_3<Kernel>` and `Sphere_3<Kernel>`
- `Line_3<Kernel>` and `Iso_cuboid_3<Kernel>`
- `Ray_3<Kernel>` and `Iso_cuboid_3<Kernel>`
- `Segment_3<Kernel>` and `Iso_cuboid_3<Kernel>`
- `Iso_cuboid_3<Kernel>` and `Iso_cuboid_3<Kernel>`.
*/
bool do_intersect(Type1<Kernel> obj1, Type2<Kernel> obj2);
/// @}



/*!
\addtogroup intersection_grp

\brief
\details Depending on which \cgal kernel is used, different overloads of this global
function are available.

\cgalHeading{Notes on Backward Compatibility}

The \ref intersection_grp function used to return an `Object`, but starting with
\cgal 4.2 the
return type is determined by a metafunction defined by the kernel. To
preserve backward compatibility `Object` can be
constructed from the new return types implicitly, but switching to the
new style is recommended. To enable the old style without any overhead,
the macro \link CGAL_INTERSECTION_VERSION `CGAL_INTERSECTION_VERSION` \endlink must be defined to
`1` before any \cgal header is included.

\sa \ref upgrading_object
\sa \ref do_intersect_grp
\sa CGAL_INTERSECTION_VERSION
*/

/*!
\addtogroup intersection_linear_grp
\ingroup intersection

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

The return type can be obtained through `CGAL::cpp11::result_of<Kernel::Intersect_2(A, B)>::%type`.
It is equivalent to `boost::optional< boost::variant< T... > >`, the last column in the table providing the template parameter pack.

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
</DIV>

\cgalHeading{3D Intersections}

The return type can be obtained through `CGAL::cpp11::result_of<Kernel::Intersect_3(A, B)>::%type`.
It is equivalent to `boost::optional< boost::variant< T... > >`, the last column in the table providing the template parameter pack.

<DIV ALIGN="CENTER">
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

\cgalHeading{Examples}

The following examples demonstrate the most common use of
`intersection()` functions with the 2D and 3D Linear %Kernel.

In the first two examples we intersect a segment and a line.
The result type can be obtained with `CGAL::cpp11::result_of`. It looks simpler
if you use a C++ compiler which supports `auto`,
but you must anyways know that the result type is a `boost::optional<boost::variant<..> >`, in order to unpack the point or segment.

<A HREF="http://www.boost.org/libs/optional/">`boost::optional`</A> comes in
as there might be no intersection. <A HREF="http://www.boost.org/libs/variant/">`boost::variant`</A> comes in
as, if there is an intersection, it is either a point or a segment.

As explained in the boost manual pages for <A HREF="http://www.boost.org/libs/variant/">`boost::variant`</A>, there are two ways to access the variants. The first examples uses `boost::get`.

\cgalExample{Kernel_23/intersection_get.cpp}

The second example uses `boost::apply_visitor`.

\cgalExample{Kernel_23/intersection_visitor.cpp}

A third example shows the use of the intersection function as a
plain function call and with `Dispatch_output_iterator`, combined with
a standard library algorithm.

\cgalExample{Kernel_23/intersections.cpp}

*/
template <typename Kernel>
cpp11::result_of<Kernel::Intersect_23(Type1, Type2)>::type
intersection(Type1<Kernel> obj1, Type2<Kernel> obj2);

/*!
returns the intersection of 3 planes, which can be a
point, a line, a plane, or empty.
*/
template <typename Kernel>
boost::optional< boost::variant< Point_3, Line_3, Plane_3 > >
intersection(const Plane_3<Kernel>& pl1,
             const Plane_3<Kernel>& pl2,
             const Plane_3<Kernel>& pl3);

/// @}
}
