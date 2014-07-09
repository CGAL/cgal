
/*!
\ingroup PkgReconstructionSimplification2Concepts
\cgalConcept

The OutputModel 
\cgalHasModel `CGAL::List_output<Tr>` 


*/

class OutputModel {
public:

/// \name Types 
/// @{


*/ 
typedef unspecified_type Output_Vertex_Iterator; 

typedef unspecified_type Output_Edge_Iterator; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns an Output_Vertex_Iterator pointing to the first vertex.
*/ 
Output_Vertex_Iterator vertices_start();

/*!
Returns an Output_Vertex_Iterator pointing beyond the last vertex.
*/ 
Output_Vertex_Iterator vertices_beyond();

/*!
Returns an Output_Edge_Iterator pointing to the first edge.
*/ 
Output_Edge_Iterator edges_start();

/*!
Returns an Output_Edge_Iterator pointing beyond the last edge.
*/ 
Output_Edge_Iterator edges_beyond();


/// @}

}; /* end OutputModel */

