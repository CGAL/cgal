namespace CGAL {

/*!
\ingroup PkgHDS_HDS

The class `HalfedgeDS_list` is a model for the `HalfedgeDS` concept. 
`HalfedgeDS_list` is a list-based representation with bidirectional 
iterators that supports removal. 

\cgalModels `HalfedgeDS<Traits,Items,Alloc>` 

\sa `CGAL::HalfedgeDS_default` 
\sa `CGAL::HalfedgeDS_vector` 
\sa `HalfedgeDSItems` 
\sa `CGAL::Polyhedron_3<Traits>` 
\sa `CGAL::HalfedgeDS_items_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_const_decorator<HDS>` 

\cgalHeading{Implementation}

`HalfedgeDS_list` uses internally the `CGAL::In_place_list` container class. 
The copy constructor and the assignment operator need \f$ O(n)\f$ time with 
\f$ n\f$ the total number of vertices, halfedges, and faces. 

`CGAL_ALLOCATOR(int)` is used as default argument for the 
`Alloc` template parameter. 

*/
template< typename Traits, typename HalfedgeDSItems, typename Alloc >
class HalfedgeDS_list {
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

/*! \name Operations 
  Besides operations required from the concept 
  `HalfedgeDS<Traits,Items,Alloc>`, this class supports 
  additionally: 
*/
/// @{

/*!
inserts elements in the range [`first, last`) before position 
`target` and removes the elements from `source`. It takes 
constant time if `&source == &``hds`; otherwise, it takes linear 
time in the size of the range. \pre [`first, last`) is a valid range in `source`. `target` is not in the range [`first, last`). 
*/ 
void vertices_splice( Vertex_iterator target, Self &source, 
Vertex_iterator first, Vertex_iterator last); 

/*!
inserts elements in the range [`first, last`) before position 
`target` and removes the elements from `source`. It takes 
constant time if `&source == &``hds`; otherwise, it takes linear 
time in the size of the range. \pre [`first, last`) is a valid range in `source`. `target` is not in the range [`first, last`). 
*/ 
void halfedges_splice( Halfedge_iterator target, Self &source, 
Halfedge_iterator first, Halfedge_iterator last); 

/*!
inserts elements in the range [`first, last`) before position 
`target` and removes the elements from `source`. It takes 
constant time if `&source == &``hds`; otherwise, it takes linear 
time in the size of the range. \pre [`first, last`) is a valid range in `source`. `target` is not in the range [`first, last`). 
*/ 
void faces_splice( Face_iterator target, Self &source, 
Face_iterator first, Face_iterator last); 

/// @}

}; /* end HalfedgeDS_list */
} /* end namespace CGAL */
