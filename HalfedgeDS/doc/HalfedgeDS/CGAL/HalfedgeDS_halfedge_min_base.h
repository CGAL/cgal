namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_VHF

The class `HalfedgeDS_halfedge_min_base` is a model of the `HalfedgeDSHalfedge`
concept. `Refs` is an instantiation of a `HalfedgeDS`.
It is equivalent to `CGAL::HalfedgeDS_halfedge_base< Refs,
CGAL::Tag_false, CGAL::Tag_false, CGAL::Tag_false>`.
The class contains support for the next and the opposite pointer and
the required type definitions. It can be used for deriving own halfedges.

\cgalModels `HalfedgeDSHalfedge`

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_min_items`
\sa `CGAL::HalfedgeDS_vertex_min_base<Refs>`
\sa `CGAL::HalfedgeDS_face_min_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>`

*/
template< typename Refs >
class HalfedgeDS_halfedge_min_base {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
HalfedgeDS_halfedge_min_base();

/// @}

}; /* end HalfedgeDS_halfedge_min_base */
} /* end namespace CGAL */
