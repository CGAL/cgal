
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Delaunay_triangulation_cell_base_3` is a model of the concept 
`DelaunayTriangulationCellBase_3`. It is the default cell base class 
of Delaunay triangulations. 


\tparam Traits must be a model of `DelaunayTriangulationTraits_3`. 

\tparam Cb must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by 
`Triangulation_cell_base_3<TriangulationTraits_3>`. 

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `DelaunayTriangulationCellBase_3` 
\sa `DelaunayTriangulationTraits_3` 
\sa `CGAL::Delaunay_triangulation_3` 
\sa `CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3`

*/

template< typename TriangulationTraits_3, typename Cb >
class Delaunay_triangulation_cell_base_3 : public Cb {
public:

/// \name Types 
/// @{
typedef TriangulationTraits_3::Point_3 Point_3;
/// @}

/*! \name Access function 

As a model of the concept `DelaunayTriangulationCellBase_3`, 
`Delaunay_triangulation_cell_base_3` 
provides a `circumcenter()` member fonction. 
*/

/// @{
/*! 
Returns the circumcenter of the cell.
*/ 
const Point_3& circumcenter(
  const TriangulationTraits_3& gt = TriangulationTraits_3()) const;

/// @}

}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
