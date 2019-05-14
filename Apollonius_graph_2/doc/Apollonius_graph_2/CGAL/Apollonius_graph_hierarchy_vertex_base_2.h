
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2

The class `Apollonius_graph_hierarchy_vertex_base_2` provides a model for the 
`ApolloniusGraphHierarchyVertexBase_2` concept, which is the 
vertex base required by the 
`Apollonius_graph_hierarchy_2<Gt,Agds>` class. The class 
`Apollonius_graph_hierarchy_vertex_base_2` is templated by a class `Agvb` which must be a model 
of the `ApolloniusGraphVertexBase_2` concept. 

\cgalModels `ApolloniusGraphHierarchyVertexBase_2`

\sa `ApolloniusGraphVertexBase_2` 
\sa `ApolloniusGraphHierarchyVertexBase_2` 
\sa `CGAL::Apollonius_graph_vertex_base_2<Gt,StoreHidden>` 

*/
template< typename Agvb >
class Apollonius_graph_hierarchy_vertex_base_2 : Agvb {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Apollonius_graph_hierarchy_vertex_base_2(); 

/*!
Constructs a vertex associated with the site `s` and 
embedded at the center of `s`. 
*/ 
Apollonius_graph_hierarchy_vertex_base_2(Site_2 s); 

/*!
Constructs a vertex associated with 
the site `s`, embedded at the center of `s`, 
and pointing to the face associated with the face 
handle `f`. 
*/ 
Apollonius_graph_vertex_base_2(Site_2 s, Face_handle f); 

/// @}

}; /* end Apollonius_graph_hierarchy_vertex_base_2 */
} /* end namespace CGAL */
