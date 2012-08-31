
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexFaceClasses

\deprecated This class is deprecated since \cgal 3.6. Its functionality is now 
transparently added when using the `Fast_location` tag as the 
`LocationPolicy` template parameter in `Delaunay_triangulation_3`, 
instead of `Triangulation_hierarchy_3`. 

This class is designed to be used as the vertex base class for 
`Triangulation_hierarchy_3`. 

It inherits from its parameter `TriangulationVertexBase_3`, and adds the 
requirements in order to match the concept 
`TriangulationHierarchyVertexBase_3`, it does so by storing two 
`Vertex_handle`s. This design allows to use either a vertex base class 
provided by \cgal, or a user customized vertex base with additional 
functionalities. 

Parameters 
-------------- 

It is parameterized by a model of the concept `TriangulationVertexBase_3`. 

\models ::TriangulationHierarchyVertexBase_3 

\sa `CGAL::Triangulation_hierarchy_3` 
\sa `CGAL::Triangulation_vertex_base_3` 
\sa `CGAL::Triangulation_vertex_base_with_info_3` 

*/
template< typename TriangulationVertexBase_3 >
class Triangulation_hierarchy_vertex_base_3 : public TriangulationVertexBase_3 {
public:

/// @}

}; /* end Triangulation_hierarchy_vertex_base_3 */
} /* end namespace CGAL */
