namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3FunctionsIO

writes a manifold or non-manifold surface
reconstructed by `make_surface_mesh()` in the OFF file format.

In case the surface is manifold the triangles can be oriented.


\tparam SurfaceMeshComplex_2InTriangulation_3 must be a model of the `SurfaceMeshComplex_2InTriangulation_3` concept.

\returns `true` if the surface is manifold and orientable.


\param os stream in which to write.
\param c2t3 Input surface.
\param options an int that is the binary union of values of `Surface_mesher::IO_option`.

\returns `true` if the surface could be written to the stream.

\see \ref IOStreamOFF
\sa `CGAL::output_surface_facets_to_polyhedron()`
*/

template <class SurfaceMeshComplex_2InTriangulation_3>
bool output_surface_facets_to_off (std::ostream& os,
                                   const  SurfaceMeshComplex_2InTriangulation_3& c2t3,
                                   int options =
                                   Surface_mesher::IO_ORIENT_SURFACE);

 namespace Surface_mesher {
/*!
\ingroup PkgSurfaceMesher3FunctionsIO
 Options for function `output_surface_facets_to_off()`.
*/
    enum IO_option { NO_OPTION = 0,
                     IO_ORIENT_SURFACE = 1,
                     IO_VERBOSE = 2 };

  } /* namespace Surface_mesher */
} /* namespace CGAL */

