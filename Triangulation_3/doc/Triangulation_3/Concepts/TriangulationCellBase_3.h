
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The cell base required by the basic triangulation does not need to store any 
geometric information, so only the requirements of the triangulation data 
structure apply. 

However, for the Delaunay and Regular triangulations, the ability to store the circumcenter 
is provided (for optimization reasons), hence an additional requirement only 
in this case, and only when the dual functions are called. 

\cgalRefines `TriangulationDSCellBase_3`

\cgalHasModel CGAL::Triangulation_cell_base_3 
\cgalHasModel CGAL::Triangulation_cell_base_with_info_3 
\cgalHasModel CGAL::Triangulation_cell_base_with_circumcenter_3 

\sa `TriangulationVertexBase_3` 

*/

class TriangulationCellBase_3 {
public:
	
/// \name Types 
/// @{

/*!
Must be the same as the point type `TriangulationTraits_3::Point_3` 
defined by the geometric traits class of the triangulation. 
*/ 
typedef unspecified_type Point_3; 
/// @} 

/// \name Access functions
/// @{
/*!
Returns the circumcenter. 
*/ 
const Point_3& circumcenter( 
const TriangulationTraits_3&gt = TriangulationTraits_3()) const; 
/// @} 

}; /* end TriangulationCellBase_3 */

