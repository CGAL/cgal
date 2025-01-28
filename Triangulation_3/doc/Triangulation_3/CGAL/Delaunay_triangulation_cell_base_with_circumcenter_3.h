
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

\tparam Traits is the geometric traits class and must be a model of `DelaunayTriangulationTraits_3`.

\tparam Cb is a cell base class from which
`Delaunay_triangulation_cell_base_with_circumcenter_3` derives. Cb should
be a model of `DelaunayTriangulationCellBase_3`.
It has the default value `Delaunay_triangulation_cell_base_3<Traits>`.

\cgalModels{DelaunayTriangulationCellBase_3}

\sa `CGAL::Delaunay_triangulation_cell_base_3`

*/
template< typename Traits, typename Cb >
class Delaunay_triangulation_cell_base_with_circumcenter_3 : public Cb {
public:

/// \name Types
/// @{
typedef Traits::Point_3 Point;
/// @}

/*! \name Access function

As a model of the concept `DelaunayTriangulationCellBase_3`,
`Delaunay_triangulation_cell_base_3`
provides a `circumcenter()` member function.

If it has already been computed in the past, the cached value is returned.
*/

/// @{

/*!
Computes the circumcenter of the tetrahedron, or retrieves it if already
computed
*/
const Point& circumcenter(Traits&gt = Traits()) const;

/// @}

}; /* end Delaunay_triangulation_cell_base_with_circumcenter_3 */
} /* end namespace CGAL */
