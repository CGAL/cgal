
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_vertex_base_3` is a model of the concept
`SurfaceMeshVertexBase_3`.
It is designed to serve as vertex base class
in a triangulation class `Tr`
plugged in a `Surface_mesh_complex_2_in_triangulation_3<Tr>`
class.

\tparam Gt is the geometric traits class.

\tparam Vb must be a model of the concept `TriangulationVertexBase_3`
and
defaults to `Triangulation_vertex_base_3 <Gt>`.

\cgalModels `SurfaceMeshVertexBase_3`

\sa `SurfaceMeshComplex_2InTriangulation_3`
\sa `Surface_mesh_complex_2_in_triangulation_3<Tr>`
\sa `SurfaceMeshTriangulation_3`
\sa `make_surface_mesh`

*/
template< typename Gt, typename Vb >
class Surface_mesh_vertex_base_3 : public Vb {
public:

}; /* end Surface_mesh_vertex_base_3 */
} /* end namespace CGAL */
