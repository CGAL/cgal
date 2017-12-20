
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

\tparam Traits is the geometric traits class and must be a model of `TriangulationTraits_3`.

\tparam Cb is a cell base class from which 
`Triangulation_cell_base_with_circumcenter_3` derives. Cb should
be a model of  `TriangulationCellBase_3` or `DelaunayTriangulationCellBase_3`. 
It has the default value `Triangulation_cell_base_3<Traits>`.

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 

*/
template< typename Traits, typename Cb >
class Triangulation_cell_base_with_circumcenter_3 : public Cb {
public:
	
/// \name Types 
/// @{
typedef Traits::Point_3 Point;
/// @}

/// \name Access Functions 
/// @{

/*!
Computes the circumcenter of the tetrahedron, or retrieves it if already 
computed
*/ 
const Point& circumcenter(const Traits& gt = Traits()) const;

/// @}

}; /* end Triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
