namespace CGAL {

/*!
\ingroup PkgKernelDFunctions

checks whether `obj1` and `obj2` intersect. Two objects `obj1` and
`obj2` intersect if there is a point `p` that is part of both `obj1`
and `obj2`.  The intersection region of those two objects is defined
as the set of all points `p` that are part of both `obj1` and
`obj2`.

\pre The objects are of the same dimension.

The types `Type1` and `Type2` can be any of the following:

- `Point_d<R>`
- `Line_d<R>`
- `Ray_d<R>`
- `Segment_d<R>`
- `Hyperplane_d<R>`

\sa `intersection`
*/
bool do_intersect(Type1<R> obj1, Type2<R> obj2);

/*!
\ingroup PkgKernelDFunctions

returns the intersection between `f1` and `f2`.

\pre The objects are of the same dimension.

The same functionality is also available through the functor `Kernel::Intersect_d`.

The following table gives the possible values for `Type1` and `Type2`.

The return type of intersecting two objects of the types `Type1` and `Type2` can be
specified through the placeholder type specifier `auto`. It is equivalent to
`std::optional< std::variant< T... > >`, the last column in the table providing
the template parameter pack.


<DIV ALIGN="CENTER">
<TABLE CELLPADDING=3 BORDER="1">
<TR>
<TH> Type1 </TH>
<TH> Type2 </TH>
<TH> `T...` </TH>
</TR>
<TR>
<TD VALIGN="CENTER" > Line_d </TD>
<TD VALIGN="CENTER" > Line_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Line_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Segment_d </TD>
<TD VALIGN="CENTER" > Line_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Segment_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Segment_d </TD>
<TD VALIGN="CENTER" > Segment_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Segment_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Ray_d </TD>
<TD VALIGN="CENTER" > Line_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Ray_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Ray_d </TD>
<TD VALIGN="CENTER" > Segment_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Segment_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Ray_d </TD>
<TD VALIGN="CENTER" > Ray_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Segment_d</TD></TR>
<TR><TD>Ray_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Hyperplane_d </TD>
<TD VALIGN="CENTER" > Line_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Line_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Hyperplane_d </TD>
<TD VALIGN="CENTER" > Ray_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Ray_d</TD></TR>
</TABLE></TD>
</TR>
<TR>
<TD VALIGN="CENTER" > Hyperplane_d </TD>
<TD VALIGN="CENTER" > Segment_d </TD>
<TD><TABLE BORDER="0">
<TR><TD>Point_d</TD></TR>
<TR><TD>Segment_d</TD></TR>
</TABLE></TD>
</TR>
</TABLE>
</DIV>

\cgalHeading{Example}

The following example demonstrates the most common use of
`intersection` routines.

\code
#include <CGAL/intersections_d.h>

template<typename R>
struct Intersection_visitor {
  typedef void result_type;
  void operator()(const Point_d<R>& p) const {
  // handle point
  }
  void operator()(const Segment_d<R>& s) const {
  // handle segment
  }
};

template <class R>
void foo(Segment_d<R> seg, Line_d<R> lin)
{
  // use auto
  auto result = intersection(seg, lin);
  if(result) { std::visit(Intersection_visitor<R>(), *result); }
  else { // no intersection
  }
}
\endcode

\sa `do_intersect`
\sa `Kernel_d::Intersect_d`
\sa <a HREF="https://www.boost.org/doc/libs/release/libs/optional/index.html">`std::optional`</a>
\sa <a HREF="https://www.boost.org/doc/html/variant.html">`std::variant`</a>

*/
decltype(auto)
intersection(Type1<R> f1, Type2<R> f2);

} /* namespace CGAL */

