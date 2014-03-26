
/*!
\ingroup PkgApolloniusGraph2Concepts
\cgalConcept

The vertex of an Apollonius graph 
included in an Apollonius graph hierarchy has to provide 
some pointers to the corresponding vertices in the 
graphs of the next and preceding levels. 
Therefore, the concept `ApolloniusGraphHierarchyVertexBase_2` 
refines the concept `ApolloniusGraphVertexBase_2`, by 
adding two vertex handles to the corresponding vertices for the 
next and previous level graphs. 

\cgalRefines `ApolloniusGraphVertexBase_2` 

\cgalHeading{Types}

`ApolloniusGraphHierarchyVertexBase_2` does not introduce any 
types in addition to those of `ApolloniusGraphVertexBase_2`. 

\cgalHasModel CGAL::Apollonius_graph_hierarchy_vertex_base_2<CGAL::Apollonius_graph_vertex_base_2<Gt,StoreHidden> > 

\sa `ApolloniusGraphDataStructure_2` 
\sa `ApolloniusGraphVertexBase_2` 
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>` 
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>` 
\sa `CGAL::Apollonius_graph_vertex_base_2<Gt,StoreHidden>` 
\sa `CGAL::Apollonius_graph_hierarchy_vertex_base_2<Agvb>` 

*/

class ApolloniusGraphHierarchyVertexBase_2 {
public:

/// \name Creation 
/// @{

/*!
Default 
constructor. 
*/ 
ApolloniusGraphHierarchyVertexBase_2(); 

/*!
Constructs a vertex associated with the site `s` and 
embedded at the center of `s`. 
*/ 
ApolloniusGraphHierarchyVertexBase_2(Site_2 s).; 

/*!
Constructs a vertex associated with 
the site `s`, embedded at the center of `s`, 
and pointing to face `f`. 
*/ 
ApolloniusGraphHierarchyVertexBase_2(Site_2 s, Face_handle f).; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns a handle to the corresponding 
vertex of the next level Apollonius graph. If such a vertex does not 
exist `Vertex_handle(NULL)` is returned. 
*/ 
Vertex_handle up(); 

/*!
Returns a handle to the corresponding 
vertex of the previous level Apollonius graph. 
*/ 
Vertex_handle down(); 

/*!
Sets the handle for the 
vertex of the next level Apollonius graph. 
*/ 
void set_up(Vertex_handle u); 

/*!
Sets the handle for the 
vertex of the previous level Apollonius graph; 
*/ 
void set_down(Vertex_handle d); 

/// @}

}; /* end ApolloniusGraphHierarchyVertexBase_2 */

