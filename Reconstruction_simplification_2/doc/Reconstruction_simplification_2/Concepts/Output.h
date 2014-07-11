
/*!
\ingroup PkgReconstructionSimplification2Concepts
\cgalConcept

The OutputModule provides a Concept which allows the user to access the simplified 
shape in a versatile way. 



\cgalHasModel `CGAL::List_output<Tr>` 
\cgalHasModel `CGAL::Off_output<Tr>` 
\cgalHasModel `CGAL::Tds_output<Tr>` 


*/

class OutputModule {
public:

/// \name Types 
/// @{


/*!
Output_Iterator for accessing the isolated vertices.
*/ 
typedef unspecified_type Output_Vertex_Iterator; 

/*!
Output_Iterator for accessing the reconstructed edges.
*/ 
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



}; /* end OutputModule */

