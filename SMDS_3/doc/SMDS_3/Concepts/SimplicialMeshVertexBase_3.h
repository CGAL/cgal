/*!
\ingroup PkgSMDS3Concepts
\cgalConcept

The concept `SimplicialMeshVertexBase_3` describes the requirements
for the `Vertex` type of the triangulation
used in the 3D simplicial mesh data structure. The type `SimplicialMeshVertexBase_3`
refines both
the concept `RegularTriangulationVertexBase_3` and the concept `SurfaceMeshVertexBase_3`.
It provides additional members to store and retrieve
information about the location of the vertex with respect
to the input domain describing the domain to be discretized.
More specifically, the concept `SimplicialMeshVertexBase_3` provides read-write access
to an integer representing the dimension of the lowest dimensional face
of the input 3D complex on which the vertex lies,
and to an index characteristic of this face.

\cgalRefines `RegularTriangulationVertexBase_3`
\cgalRefines `SurfaceMeshVertexBase_3`

\cgalHasModel `CGAL::Mesh_vertex_base_3<Gt,MD,Vb>`
\cgalHasModel `CGAL::Simplicial_mesh_vertex_base_3<Gt>`

*/

class MeshVertexBase_3 {
public:

/// \todo update this file from here

/// \name Types
/// @{

/*!
Index type. Must match the type `MeshDomain_3::Index`.
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

/// @}

/*! \name Internal
These functions are used internally. The user is
not encouraged to use them directly as they may change in the future.
*/
/// @{

/*!
Returns a boolean, used for feature edges protection.
*/
bool is_special();

/*!
Sets the special aspect of the vertex.
*/
void set_special(bool);

/*!

*/
FT meshing_info() const;

/*!

*/
void set_meshing_info(FT);

/*!

*/
Vertex_handle next_intrusive() const;

/*!

*/
void set_next_intrusive(Vertex_handle);

/*!

*/
Vertex_handle previous_intrusive() const;

/*!

*/
void set_previous_intrusive(Vertex_handle);

/// Get the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
unsigned int erase_counter() const;

/// Sets the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
void set_erase_counter(unsigned int c);

/// Increments the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
void increment_erase_counter();
/// @}

}; /* end MeshVertexBase_3 */
