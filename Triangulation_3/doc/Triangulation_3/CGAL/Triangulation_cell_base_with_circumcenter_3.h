
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_cell_base_with_circumcenter_3` is a model of the 
concept `TriangulationCellBase_3`, the base cell of a 3D-triangulation. 
It provides an easy way to cache the computation of the circumcenter of 
tetrahedra. 
Note that input/output operators discard this additional information. 

All functions modifying the vertices of the cell, invalidate the cached 
circumcenter. 

\tparam TriangulationTraits_3 is the geometric traits class. 

\tparam TriangulationCellBase_3 is a cell base class from which 
`Triangulation_cell_base_with_circumcenter_3` derives. 
It has the default value `Triangulation_cell_base_3<TriangulationTraits_3>`. 

\cgalModels `TriangulationCellBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 

*/
template< typename TriangulationTraits_3, typename TriangulationCellBase_3 >
class Triangulation_cell_base_with_circumcenter_3 : public TriangulationCellBase_3 {
public:
	
/// \name Types 
/// @{
typedef TriangulationTraits_3::Point_3 Point_3;
/// @}

/// \name Access Functions 
/// @{

/*!
Computes the circumcenter of the tetrahedron, or retrieve it if already 
computed. 
*/ 
const Point_3& circumcenter( 
const TriangulationTraits_3&gt = TriangulationTraits_3()) const; 

/// @}

}; /* end Triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
