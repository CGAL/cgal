
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Delaunay_triangulation_cell_base_3` is a model of the concept 
`DelaunayTriangulationCellBase_3`. It is the default cell base class 
of Delaunay triangulations. 

\tparam DelaunayTriangulationTraits_3 is the geometric traits class.

\tparam Cb must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by 
`Triangulation_cell_base_3<DelaunayTriangulationTraits_3>`. 

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `DelaunayTriangulationCellBase_3` 
\sa `CGAL::Delaunay_triangulation_3` 
\sa `CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3`

*/

template< typename DelaunayTriangulationTraits_3, typename Cb >
class Delaunay_triangulation_cell_base_3 : public Cb {
public:

/// \name Types 
/// @{
typedef DelaunayTriangulationTraits_3::Point_3 Point_3;
/// @}

/*! \name Access function 

As a model of the concept `DelaunayTriangulationCellBase_3`, 
`Delaunay_triangulation_cell_base_3` 
provides a `circumcenter()` member fonction. 
*/

/// @{
/*! 
Returns the circumcenter of the cell
*/ 
const Point_3& circumcenter(
  const DelaunayTriangulationTraits_3& gt = DelaunayTriangulationTraits_3()) const;

/// @}

}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
