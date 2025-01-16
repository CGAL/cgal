
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `SurfaceMeshVertexBase_3` describes the vertex base type
of the three dimensional triangulation used
to embed the surface mesh.

More precisely,
the first template parameter `SurfaceMeshC2T3` of the function template
`CGAL::make_surface_mesh()`
is a model of the concept
`SurfaceMeshComplex_2InTriangulation_3`
which describes a data structure to store
a pure two dimensional complex
embedded in a three dimensional triangulation.
In particular, the type `SurfaceMeshC2T3` is required to provide
a three dimensional triangulation type
`SurfaceMeshC2T3::Triangulation_3`
The concept `SurfaceMeshVertexBase_3` describes the vertex base type
required in this triangulation type.

\cgalRefines{TriangulationVertexBase_3}
The surface mesher algorithm issues frequent queries about
the status of the vertices with respect to the two dimensional complex
that represents the current surface approximation.  The class
SurfaceMeshVertexBase_3 offers a caching mechanism to answer more
efficiently these queries. The caching mechanism includes two cached
integers, which, when they are valid, store respectively the number of
complex facets incident to the vertex and the number of connected
components of the adjacency graph of those facets.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Surface_mesh_vertex_base_3<Gt,Vb>}
\cgalHasModels{CGAL::Surface_mesh_default_triangulation_3::Vertex}
\cgalHasModelsEnd

\sa `SurfaceMeshComplex_2InTriangulation_3`
\sa `CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>`
\sa `CGAL::Surface_mesh_default_triangulation_3`

*/

class SurfaceMeshVertexBase_3 {
public:

/// \name Operations
/// @{

/*!
Returns `true` if the cache is valid.
*/
bool is_c2t3_cache_valid();

/*!
Invalidates the cache.
*/
void invalidate_c2t3_cache();

/*!
Returns the cached number of facets of the complex incident to the vertex.
*/
int cached_number_of_incident_facets();

/*!
This method concerns the adjacency
graph of the facets of the complex incident to the vertex
and returns a cached value for the number of connected components this graph.
*/
int cached_number_of_components();

/// @}

}; /* end SurfaceMeshVertexBase_3 */

