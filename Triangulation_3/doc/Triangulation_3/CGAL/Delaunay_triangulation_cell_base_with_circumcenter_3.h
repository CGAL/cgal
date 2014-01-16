
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Delaunay_triangulation_cell_base_with_circumcenter_3` derives from 
`Cb`, a cell base class of a 3D triangulation.
It provides an easy way to cache the computation of the circumcenters of 
tetrahedra. 
Note that input/output operators discard this additional information. 

All functions modifying the vertices of the cell invalidate the cached 
circumcenter. 

\tparam TriangulationTraits_3 is the geometric traits class. It should be a model
of `DelaunayTriangulationTraits_3`.

\tparam Cb is a cell base class from which 
`Delaunay_triangulation_cell_base_with_circumcenter_3` derives. Cb should
be a model of `DelaunayTriangulationCellBase_3`. 
It has the default value `Triangulation_cell_base_3<DelaunayTriangulationTraits_3>`. 

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 
\sa `CGAL::Delaunay_triangulation_cell_base_3`

*/
template< typename TriangulationTraits_3, typename Cb >
class Delaunay_triangulation_cell_base_with_circumcenter_3 : public Cb {
public:
	
/// \name Types 
/// @{
typedef TriangulationTraits_3::Point_3 Point_3;
/// @}

/// \name Access Functions 
/// @{

/*!
Computes the circumcenter of the tetrahedron, or retrieves it if already 
computed. 
*/ 
const Point_3& circumcenter( 
const TriangulationTraits_3&gt = TriangulationTraits_3()) const; 

/// @}

}; /* end Delaunay_triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
