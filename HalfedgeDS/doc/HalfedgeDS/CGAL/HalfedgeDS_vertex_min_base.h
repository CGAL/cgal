namespace CGAL {

/*!
\ingroup PkgHDS_VHF

The class `HalfedgeDS_vertex_min_base` is a model of the `HalfedgeDSVertex` 
concept. `Refs` is an instantiation of a `HalfedgeDS`. It is 
equivalent to `HalfedgeDS_vertex_base< Refs, CGAL::Tag_false>`. 
It is empty besides the required type definitions. It can be used for deriving 
own vertices. 

\cgalModels `HalfedgeDSVertex`

\sa `HalfedgeDS<Traits,Items,Alloc>` 
\sa `HalfedgeDSItems` 
\sa `PolyhedronItems_3` 
\sa `CGAL::HalfedgeDS_min_items` 
\sa `CGAL::HalfedgeDS_halfedge_min_base<Refs>` 
\sa `CGAL::HalfedgeDS_face_min_base<Refs>` 
\sa `CGAL::HalfedgeDS_vertex_base<Refs>` 

*/
template< typename Refs >
class HalfedgeDS_vertex_min_base {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
HalfedgeDS_vertex_min_base(); 

/// @}

}; /* end HalfedgeDS_vertex_min_base */
} /* end namespace CGAL */
