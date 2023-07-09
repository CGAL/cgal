namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_VHF

The class `HalfedgeDS_halfedge_base` is a model of the `HalfedgeDSHalfedge`
concept. `Refs` is an instantiation of a `HalfedgeDS`.
The full declaration states four template parameters:

<TABLE border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<span class="mbox"></span>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`template <`
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`class Refs,`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`class Tag_prev `
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`= CGAL::Tag_true,`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`class Tag_vertex `
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`= CGAL::Tag_true,`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`class Tag_face `
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`= CGAL::Tag_true>`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`class HalfedgeDS_halfedge_base;`

</TABLE>

If `Tag_prev` \f$ \equiv\f$ `CGAL::Tag_true` a reference to the previous
halfedge is supported.

If `Tag_vertex` \f$ \equiv\f$ `CGAL::Tag_true` an incident vertex is
supported.

If `Tag_face` \f$ \equiv\f$ `CGAL::Tag_true` an incident face is
supported.

In all cases, a reference to the next halfedge and to the opposite
halfedge is supported.

\cgalModels{HalfedgeDSHalfedge}

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::HalfedgeDS_vertex_base<Refs>`
\sa `CGAL::HalfedgeDS_face_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_min_base<Refs>`

*/
template< typename Refs >
class HalfedgeDS_halfedge_base {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
HalfedgeDS_halfedge_base();

/// @}

}; /* end HalfedgeDS_halfedge_base */
} /* end namespace CGAL */
