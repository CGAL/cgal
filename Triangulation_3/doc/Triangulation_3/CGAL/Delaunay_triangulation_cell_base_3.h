
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Delaunay_triangulation_cell_base_3` is a model of the concept
`DelaunayTriangulationCellBase_3`. It is the default cell base class
of Delaunay triangulations.

\tparam Traits is the geometric traits class and must be a model of `DelaunayTriangulationTraits_3`.

\tparam Cb must be a model of `TriangulationCellBase_3`.
By default, this parameter is instantiated by
`Triangulation_cell_base_3<Traits>`.

\cgalModels `DelaunayTriangulationCellBase_3`

\sa `DelaunayTriangulationCellBase_3`
\sa `CGAL::Delaunay_triangulation_3`
\sa `CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3`

*/

template< typename Traits, typename Cb >
class Delaunay_triangulation_cell_base_3 : public Cb {
public:

/// \name Types
/// @{
typedef Traits::Point_3 Point;
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
const Point& circumcenter(const Traits& gt = Traits()) const;

/// @}

}; /* end Delaunay_triangulation_cell_base_3 */
} /* end namespace CGAL */
