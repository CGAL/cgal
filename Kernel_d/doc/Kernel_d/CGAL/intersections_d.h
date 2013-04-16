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

returns the intersection result of `f1` and `f2` by means of
the polymorphic wrapper type `Object`. The returned object can be
tested for the intersection result and assigned by means of the
`object_cast` function. 

\pre The objects are of the same dimension.

The possible value for types `Type1` and `Type2` and 
the possible return values wrapped in `Object` are the following: 

<DIV ALIGN="CENTER"> 
<TABLE CELLPADDING=3 BORDER="1"> 
<TR> <TH> Type1 </TH> 
<TH> Type2 </TH> 
<TH> Return Type </TH> 
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

template <class R> 
void foo(Segment_d<R> seg, Line_d<R> lin) 
{ 
  Object result = intersection(seg, lin); 
  if (const Point_d<R> *ipnt = object_cast<Point_d<R> >(&result) ) { 
  // handle the point intersection case with *ipnt. 
  } else if (const Segment_d<R> *iseg = object_cast<Segment_d<R> >(&result) ) { 
  // handle the segment intersection case with *iseg. 
  } else { 
  // handle the no intersection case. 
  } 
} 
\endcode

\sa `do_intersect`, \sa `Kernel_d::Intersect_d` 

*/
Object intersection(Type1<R> f1, Type2<R> f2);

} /* namespace CGAL */

