
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The base cell of a Delaunay triangulation 
has to be a model 
of the concept `DelaunayTriangulationCellBase_3`, which refines 
the concept `TriangulationCellBase_3` by adding 
in the cell an operator that computes its circumcenter. 

For the Delaunay triangulation,
the model `CGAL::Triangulation_cell_base_with_circumcenter_3`
gives the ability to store the circumcenter (for optimization reasons),
hence this additional requirement. 


\cgalRefines `TriangulationCellBase_3`

\cgalHasModel CGAL::Delaunay_triangulation_cell_base_3 
\cgalHasModel CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3

\sa `TriangulationCellBase_3` 

*/

class DelaunayTriangulationCellBase_3 {
public:

/// \name Access functions
/// @{
/*!
Returns the circumcenter of the cell. 
`DelaunayTriangulationTraits_3` is the geometric traits class of the triangulation.
*/ 
const DelaunayTriangulationTraits_3::Point_3& circumcenter( 
const DelaunayTriangulationTraits_3&gt = DelaunayTriangulationTraits_3()) const; 
/// @} 


}; /* end DelaunayTriangulationCellBase_3 */
