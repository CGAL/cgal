
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

\deprecated This concept is deprecated since \cgal 3.6, as the only class using 
it, `Triangulation_hierarchy_3` has been deprecated as well. 

The vertex base used by `Triangulation_hierarchy_3` must provide 
access to two vertex handles for linking between the levels of the hierarchy. 

\cgalRefines `TriangulationVertexBase_3` 

\cgalHasModel CGAL::Triangulation_hierarchy_vertex_base_3 

*/

class TriangulationHierarchyVertexBase_3 {
public:

/// \name Access Functions 
/// @{

/*! 
Returns the `Vertex_handle` pointing to the level above. 
*/ 
Vertex_handle up() const; 

/*! 
Returns the `Vertex_handle` pointing to the level below. 
*/ 
Vertex_handle down() const; 

/// @} 

/// \name Setting 
/// @{

/*! 
Sets the `Vertex_handle` pointing to the level above to `v`. 
*/ 
void set_up(Vertex_handle v); 

/*! 
Sets the `Vertex_handle` pointing to the level below to `v`. 
*/ 
void set_down(Vertex_handle v); 

/// @}

}; /* end TriangulationHierarchyVertexBase_3 */

