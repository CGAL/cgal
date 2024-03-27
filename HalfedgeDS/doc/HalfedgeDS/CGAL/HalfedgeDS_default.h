namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_HDS

The class `HalfedgeDS_default` is a model for the `HalfedgeDS` concept. The
second template parameter `HalfedgeDSItems` has a default argument
`CGAL::HalfedgeDS_items_2`. The third template parameter `Alloc`
uses the \cgal default allocator as default setting. `HalfedgeDS_default` is a
list-based representation with bidirectional iterators that supports
removal.

\cgalModelsBareBegin
\cgalModelsBare{`HalfedgeDS<Traits,Items,Alloc>`}
\cgalModelsBareEnd

\sa `CGAL::HalfedgeDS_list`
\sa `CGAL::HalfedgeDS_vector`
\sa `HalfedgeDSItems`
\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::Polyhedron_3<Traits>`
\sa `CGAL::HalfedgeDS_items_decorator<HDS>`
\sa `CGAL::HalfedgeDS_decorator<HDS>`
\sa `CGAL::HalfedgeDS_const_decorator<HDS>`

\cgalHeading{Implementation}

Currently, `HalfedgeDS_default` is derived from `CGAL::HalfedgeDS_list<Traits>`.
The copy constructor and the assignment operator need \cgalBigO{n} time with
\f$ n\f$ the total number of vertices, halfedges, and faces.

*/
template< typename Traits,
          typename HalfedgeDSItems = CGAL::Halfedge_DS_items_2,
          typename Alloc = CGAL_ALLOCATOR(int) >
class HalfedgeDS_default {
public:

/// \name Types
/// @{

/*!

*/
typedef bidirectional_iterator_tag iterator_category;

/*!

*/
typedef CGAL::Tag_true Supports_removal;

/// @}

}; /* end HalfedgeDS_default */
} /* end namespace CGAL */
