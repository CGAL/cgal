
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses
\deprecated This class is deprecated since \cgal 4.4. The class
`Delaunay_triangulation_cell_base_with_circumcenter_3` should
be used instead. 
`Regular_triangulation_cell_base_with_weighted_circumcenter_3` and
`Triangulation_cell_base_with_info_3` can also be used in other frameworks.

The class `Triangulation_cell_base_with_circumcenter_3` derives from 
`Cb`, a cell base class of a 3D triangulation.
It provides an easy way to cache the computation of the circumcenters of 
tetrahedra. 
Note that input/output operators discard this additional information. 

All functions modifying the vertices of the cell, invalidate the cached 
circumcenter. 

\tparam TriangulationTraits_3 is the geometric traits class. 

\tparam Cb is a cell base class from which 
`Triangulation_cell_base_with_circumcenter_3` derives. Cb should
be a model of  `TriangulationCellBase_3` or `DelaunayTriangulationCellBase_3`. 
It has the default value `Triangulation_cell_base_3<DelaunayTriangulationTraits_3>`. 

\cgalModels `TriangulationCellBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 

*/
template< typename TriangulationTraits_3, typename Cb >
class Triangulation_cell_base_with_circumcenter_3 : public Cb {
public:
	
/// \name Types 
/// @{
typedef TriangulationTraits_3::Point_3 Point_3;
/// @}

/// \name Access Functions 
/// @{

/*!
Computes the circumcenter of the tetrahedron, or retrieves it if already 
computed
*/ 
const Point_3& circumcenter( 
const TriangulationTraits_3&gt = TriangulationTraits_3()) const; 

/// @}

}; /* end Triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
