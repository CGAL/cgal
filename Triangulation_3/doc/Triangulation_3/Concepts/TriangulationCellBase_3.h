
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The cell base required by the basic triangulation does not need to store any 
geometric information, so only the requirements of the triangulation data 
structure apply. 


\cgalRefines `TriangulationDSCellBase_3`

\cgalHasModel CGAL::Triangulation_cell_base_3 
\cgalHasModel CGAL::Triangulation_cell_base_with_info_3 

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
typedef unspecified_type Point; 

/// @} 

	
}; /* end TriangulationCellBase_3 */