
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_complex_2_in_triangulation_3` implements a data structure to store
the restricted Delaunay triangulation used by the surface mesh generator.
The restricted Delaunay triangulation is stored as a two dimensional
complex embedded in a three dimensional triangulation.

The class `Surface_mesh_complex_2_in_triangulation_3` is a model of the concept `SurfaceMeshComplex_2InTriangulation_3`
and can be plugged as the template parameter `C2T3`
in the function template `make_surface_mesh()`.

\tparam Tr has to be instantiated
with a model of the concept `SurfaceMeshTriangulation_3`.
(Any three dimensional triangulation of
\cgal is a model of `Tr`
provided that its vertex and cell base class be models
of the concept `SurfaceMeshVertexBase_3` and
`SurfaceMeshCellBase_3`
respectively.)

\cgalModels{SurfaceMeshComplex_2InTriangulation_3}

\sa `make_surface_mesh`
\sa `SurfaceMeshTriangulation_3`
\sa `CGAL::output_surface_facets_to_off()`
\sa `CGAL::output_surface_facets_to_polyhedron()`

*/
template< typename Tr, typename Edge_info = void >
class Surface_mesh_complex_2_in_triangulation_3 {
public:

}; /* end Surface_mesh_complex_2_in_triangulation_3 */
} /* end namespace CGAL */
