
/*!
\ingroup PkgReconstructionSimplification2Concepts
\cgalConcept

The OutputModule provides a Concept which allows the user to access the simplified 
shape in a versatile way. 



\cgalHasModel `CGAL::List_output<Kernel, Output_Vertex_Iterator, Output_Edge_Iterator>` 
\cgalHasModel `CGAL::Off_output<Kernel>` 
\cgalHasModel `CGAL::Tds_output<Kernel>` 


*/

class OutputModule {
public:

/*!
Extracts the solid edges and vertices from the `Reconstruction_simplification_2` module.

\param rt2 The `Reconstruction_triangulation_2` from which the solid edges and vertices are extracted.
\param nb_ignore The number of verticess to be ignored in the output.
		
*/ 
void store_marked_elements(Rt_2& rt2, int nb_ignore);


}; /* end OutputModule */

