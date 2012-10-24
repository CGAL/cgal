
/*!
\ingroup PkgTriangulation3Concepts
\cgalconcept

The cell base required by the basic triangulation does not need to store any 
geometric information, so only the requirements of the triangulation data 
structure apply. 

However, for the Delaunay triangulation, the ability to store the circumcenter 
is provided (for optimization reasons), hence an additional requirement only 
in this case, and only when the dual functions are called. 

\refines `TriangulationDSCellBase_3`

\hasModel CGAL::Triangulation_cell_base_3 
\hasModel CGAL::Triangulation_cell_base_with_info_3 
\hasModel CGAL::Triangulation_cell_base_with_circumcenter_3 

\sa `TriangulationVertexBase_3` 


*/

class TriangulationCellBase_3 {
public:
/*! 
Returns the circumcenter. 
*/ 
const DelaunayTriangulationTraits_3::Point_3& circumcenter( 
const DelaunayTriangulationTraits_3&gt = DelaunayTriangulationTraits_3()) const; 

}; /* end TriangulationCellBase_3 */

