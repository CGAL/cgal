/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalconcept

The concept `MeshVertexBase_3` describes the requirements 
for the `Vertex` type of the triangulation 
used by a 3D mesh generation process. The type `MeshVertexBase_3` refines both the concept `TriangulationVertexBase_3` 
and 
the concept `SurfaceMeshVertexBase_3`. 
It provides additional members to store and retrieve 
information about the location of the vertex with respect 
to the input domain describing the domain to be discretized. 
More specifically, the concept `MeshVertexBase_3` provides read-write access 
to an integer representing the dimension of the lowest dimensional face 
of the input 3D complex on which the vertex lies, 
and to an index characteristic of this face. 

\refines `TriangulationVertexBase_3` 
\refines `SurfaceMeshVertexBase_3`

\hasModel `CGAL::Mesh_vertex_base_3<MD,Gt,Vb>` 

\sa `CGAL::make_mesh_3` 
\sa `CGAL::refine_mesh_3` 
\sa `MeshDomain_3` 

*/

class MeshVertexBase_3 {
public:

/// \name Types 
/// @{

/*! 
Index type. Must match the type `MeshDomain_3::Index`. 
*/ 
typedef Hidden_type Index;; 

/*! 
Numerical type. 
*/ 
typedef Hidden_type FT;; 

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
These functions are used internally by mesh optimizers. The user is 
not encouraged to use them directly as they may change in the future. 
*/
/// @{

/*! 

*/ 
FT meshing_info() const; 

/*! 

*/ 
void set_meshing_info(FT); 

/// @}

}; /* end MeshVertexBase_3 */
