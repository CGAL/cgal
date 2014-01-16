
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Delaunay_triangulation_cell_base_3` is a model of the concept 
`DelaunayTriangulationCellBase_3`. It is the default cell base class 
of Delaunay triangulations. 


\tparam Traits must be a model of `DelaunayTriangulationTraits_3`. 

\tparam Cb must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by `Triangulation_cell_base_3<Traits>`. 

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `DelaunayTriangulationCellBase_3` 
\sa `DelaunayTriangulationTraits_3` 
\sa `CGAL::Delaunay_triangulation_3<Traits,Tds>` 
\sa `CGAL::Triangulation_cell_base_with_circumcenter_3<DelaunayTriangulationTraits_3, CellBase_3>`

*/

template< typename Traits, typename Cb >
class Delaunay_triangulation_cell_base_3 : public Cb {
public:

/*! \name Access function 

As a model of the concept `DelaunayTriangulationCellBase_3`, 
`Delaunay_triangulation_cell_base_3` 
provides a `circumcenter()` member fonction. 

A class for the cells of Delaunay triangulations with
cached circumcenters can simply be obtained by plugging 
`Delaunay_triangulation_cell_base_3` in the second template parameter of 
`CGAL::Triangulation_cell_base_with_circumcenter_3<DelaunayTriangulationTraits_3, CellBase_3>`.
*/

/// @{
/*! 
Returns the circumcenter of the cell.
*/ 
const Traits::Point_3& circumcenter(const Traits& gt = Traits()) const;

/// @}

}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
