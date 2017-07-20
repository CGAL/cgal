namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3FunctionsIO

 Gets reconstructed surface out of a `SurfaceMeshComplex_2InTriangulation_3` object.

 This variant exports the surface as a `FaceGraph` and appends it to `graph`.
 It must be manifold. For this purpose, you may call
 `make_surface_mesh()` with `Manifold_tag` or
 `Manifold_with_boundary_tag` parameter.

 @tparam C2T3 model of the `SurfaceMeshComplex_2InTriangulation_3` concept.
 @tparam FaceGraph a model of `MutableFaceGraph`.

 @param c2t3 an instance of a manifold `C2T3`.
 @param graph an instance of `FaceGraph`.

\sa `CGAL::output_surface_facets_to_off()`
*/
template<class C2T3, class FaceGraph>
void output_surface_facets_to_facegraph(const C2T3& c2t3, FaceGraph& graph);

} /* namespace CGAL */

