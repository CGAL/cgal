
namespace CGAL {

/*!
\ingroup PkgPolyhedron

The class `Polyhedron_min_items_3` is a minimal model of the 
`PolyhedronItems_3` concept. It provides definitions for vertices 
containing points, halfedges, and faces. The polyhedron traits class 
must provide the respective type for the point. Vertices and facets 
both do <I>not</I> contain a halfedge handle to an incident 
halfedge. 

\cgalModels `PolyhedronItems_3`

\sa `CGAL::Polyhedron_3<Traits>`
\sa `CGAL::Polyhedron_items_3` 
\sa `CGAL::HalfedgeDS_min_items` 
\sa `CGAL::HalfedgeDS_items_2` 

*/

class Polyhedron_min_items_3 {
public:

/// \name Types in Polyhedron_min_items_3::Face_wrapper<Refs,Traits>::Vertex 
/// @{

`Polyhedron_min_items_3::Vertex_wrapper<Refs,Traits>::Vertex` 
/*!

*/ 
typedef Traits::Point_3 Point; 

`Polyhedron_min_items_3::Vertex_wrapper<Refs,Traits>::Vertex` 
/*!

*/ 
typedef CGAL::Tag_true Supports_vertex_point; 

/// @} 

/// \name Types in Polyhedron_min_items_3::Face_wrapper<Refs,Traits>::Face
/// @{

/*!

*/ 
typedef CGAL::Tag_false Supports_face_plane; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Polyhedron_min_items_3(); 

/// @}

}; /* end Polyhedron_min_items_3 */
} /* end namespace CGAL */
