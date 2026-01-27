namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

This class is an alias for the current state-of-the-art Garland-Heckbert policies
(see Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).
It is currently an alias of `GarlandHeckbert_plane_and_line_policies`.

\sa `GarlandHeckbert_plane_policies`
\sa `GarlandHeckbert_probabilistic_plane_policies`
\sa `GarlandHeckbert_triangle_policies`
\sa `GarlandHeckbert_probabilistic_triangle_policies`
\sa `GarlandHeckbert_plane_and_line_policies`

*/
template<typename TriangleMesh, typename GeomTraits>
using GarlandHeckbert_policies = GarlandHeckbert_plane_and_line_policies<TriangleMesh, GeomTraits>;

} // namespace Surface_mesh_simplification
} // namespace CGAL
