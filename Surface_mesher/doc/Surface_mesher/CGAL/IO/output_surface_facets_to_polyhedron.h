namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3

`output_surface_facets_to_polyhedron()` converts a surface
reconstructed by `make_surface_mesh()` to a `Polyhedron_3<Traits>`.

Gets reconstructed surface out of a
`SurfaceMeshComplex_2InTriangulation_3` object.

This variant exports the surface as a polyhedron. It requires the
surface to be manifold. For this purpose, you may call
`make_surface_mesh`() with `Manifold_tag` or
`Manifold_with_boundary_tag` parameter.

### Template Parameters ###

`SurfaceMeshComplex_2InTriangulation_3`: model of the `SurfaceMeshComplex_2InTriangulation_3` concept. `Polyhedron`: an instance of `Polyhedron_3<Traits>`.

### Returns ###

true if the surface is manifold and orientable.

### Parameters ###

`c2t3`: Input surface. `output_polyhedron`: Output polyhedron.

\sa `CGAL::output_surface_facets_to_off`
*/
template<class SurfaceMeshComplex_2InTriangulation_3, class Polyhedron> bool output_surface_facets_to_polyhedron(const SurfaceMeshComplex_2InTriangulation_3& c2t3, Polyhedron& output_polyhedron);

} /* namespace CGAL */

