/*!
\ingroup PkgSMDS3Concepts
\cgalConcept

The concept `SimplicialMeshVertexBase_3` describes the requirements
for the `Vertex` type of the triangulation
used in a 3D simplicial mesh data structure. The type `SimplicialMeshVertexBase_3`
refines the concept `TriangulationVertexBase_3`.
It provides additional members to store and retrieve
information about the location of the vertex with respect
to the input domain describing the discretized domain.
More specifically, the concept `SimplicialMeshVertexBase_3` provides read-write access
to an integer representing the dimension of the lowest dimensional face
of the input 3D complex on which the vertex lies,
and to an index characteristic of this face.

\cgalRefines{TriangulationVertexBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_vertex_base_3}
\cgalHasModels{CGAL::Simplicial_mesh_vertex_base_3}
\cgalHasModels{CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3}
\cgalHasModelsEnd

*/

class SimplicialMeshVertexBase_3 {
public:

/// \name Types
/// @{

/*!
Index type.
*/
typedef unspecified_type Index;

/*!
Numerical type.
*/
typedef unspecified_type FT;

/// @}

/// \name Operations
/// @{

/*!
Returns the dimension of the lowest dimensional face of the input 3D complex that contains the vertex.
*/
int in_dimension() const;

/*!
Sets the dimension of the lowest dimensional face of the input 3D complex that contains the vertex.
*/
void set_dimension(int);

/*!
Returns the index of the lowest dimensional face of the input 3D complex that contains the vertex.
*/
Index index();

/*!
Sets the index of the lowest dimensional face of the input 3D complex that contains the vertex.
*/
void set_index(Index);

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


}; /* end SimplicialMeshVertexBase_3 */
