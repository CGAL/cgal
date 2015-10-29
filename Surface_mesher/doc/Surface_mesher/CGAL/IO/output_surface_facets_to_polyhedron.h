namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3FunctionsIO

converts a manifold surface reconstructed by `make_surface_mesh()` to a `Polyhedron_3<Traits>`.

The surface must be manifold. For this purpose, you may call
`make_surface_mesh()` with `Manifold_tag` or
`Manifold_with_boundary_tag` parameter.

\tparam SurfaceMeshComplex_2InTriangulation_3 must be a model of the `SurfaceMeshComplex_2InTriangulation_3` concept. 
\tparam Polyhedron must be an instance of `Polyhedron_3<Traits>`.

\returns `true` if the surface is manifold and orientable.


\param c2t3 Input surface. 
\param output_polyhedron Output polyhedron.

\sa `CGAL::output_surface_facets_to_off()`
*/
template<class SurfaceMeshComplex_2InTriangulation_3, class Polyhedron> bool output_surface_facets_to_polyhedron(const SurfaceMeshComplex_2InTriangulation_3& c2t3, Polyhedron& output_polyhedron);

} /* namespace CGAL */

