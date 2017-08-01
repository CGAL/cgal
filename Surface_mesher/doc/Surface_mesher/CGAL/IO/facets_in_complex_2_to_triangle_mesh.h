namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3FunctionsIO

 \brief Converts a manifold surface reconstructed by `make_surface_mesh()` to a `TriangleMesh`.

 This variant exports the surface as a `TriangleMesh` and appends it to `graph`.
 It must be manifold. For this purpose, you may call
 `make_surface_mesh()` with `Manifold_tag` or
 `Manifold_with_boundary_tag` parameter.

 @tparam C2T3 model of the `SurfaceMeshComplex_2InTriangulation_3` concept.
 @tparam TriangleMesh a model of `MutableFaceGraph`.

 @param c2t3 an instance of a manifold `C2T3`.
 @param graph an instance of `TriangleMesh`.

\sa `CGAL::output_surface_facets_to_off()`
*/
template<class C2T3, class TriangleMesh>
void facets_in_complex_2_to_triangle_mesh(const C2T3& c2t3, TriangleMesh& graph);

} /* namespace CGAL */

