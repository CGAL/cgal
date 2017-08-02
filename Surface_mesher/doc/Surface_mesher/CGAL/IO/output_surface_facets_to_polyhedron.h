namespace CGAL {

/*!
\deprecated Gets reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.

This variant exports the surface as a polyhedron.
It requires the surface to be manifold. For this purpose,
you may call make_surface_mesh() with Manifold_tag or Manifold_with_boundary_tag parameter.

@commentheading Template Parameters:
@param SurfaceMeshComplex_2InTriangulation_3 model of the SurfaceMeshComplex_2InTriangulation_3 concept.
@param Polyhedron an instance of CGAL::Polyhedron_3<Traits>.

@return true if the surface is manifold and orientable.
*/
template<class SurfaceMeshComplex_2InTriangulation_3, class Polyhedron> bool output_surface_facets_to_polyhedron(const SurfaceMeshComplex_2InTriangulation_3& c2t3, Polyhedron& output_polyhedron);

} /* namespace CGAL */

