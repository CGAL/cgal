
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexFaceClasses

The class `Triangulation_cell_base_with_circumcenter_3` is a model of the 
concept `TriangulationCellBase_3`, the base cell of a 3D-triangulation. 
It provides an easy way to cache the computation of the circumcenter of 
tetrahedra. 
Note that input/output operators discard this additional information. 

All functions modifying the vertices of the cell, invalidate the cached 
circumcenter. 

### Parameters ###

The first template argument is the geometric traits class 
`DelaunayTriangulationTraits_3`. 

The second template argument is a cell base class from which 
`Triangulation_cell_base_with_circumcenter_3` derives. 
It has the default value `Triangulation_cell_base_3<DelaunayTriangulationTraits_3>`. 

\models ::TriangulationCellBase_3 

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 

*/
template< typename DelaunayTriangulationTraits_3, typename TriangulationCellBase_3 >
class Triangulation_cell_base_with_circumcenter_3 : public TriangulationCellBase_3 {
public:

/// \name Access Functions 
/// @{

/*! 
Computes the circumcenter of the tetrahedron, or retrieve it if already 
computed. 
*/ 
const DelaunayTriangulationTraits_3::Point_3& circumcenter( 
const DelaunayTriangulationTraits_3&gt = DelaunayTriangulationTraits_3()) const; 

/// @}

}; /* end Triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
