
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_cell_base_3` is a model of the concept
`SurfaceMeshCellBase_3`.
It is designed to serve as cell base class
in a triangulation class `Tr`
plugged in a `Surface_mesh_complex_2_in_triangulation_3<Tr>`
class.

\tparam Gt is the geometric traits class.

\tparam Cb must be a model of the concept `DelaunayTriangulationCellBase_3`
and defaults to `Delaunay_triangulation_cell_base_3<GT>`.

\cgalModels `SurfaceMeshCellBase_3`

\sa `SurfaceMeshComplex_2InTriangulation_3`
\sa `Surface_mesh_complex_2_in_triangulation_3<Tr>`
\sa `SurfaceMeshTriangulation_3`
\sa `make_surface_mesh`
*/
template< typename Gt, typename Cb >
class Surface_mesh_cell_base_3 : public Cb {
public:

}; /* end Surface_mesh_cell_base_3 */
} /* end namespace CGAL */
