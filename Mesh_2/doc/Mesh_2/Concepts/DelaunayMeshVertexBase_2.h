


/*!
\ingroup PkgMesh2Concepts
\cgalConcept


The concept `DelaunayMeshVertexBase_2` refines the concept
`TriangulationVertexBase_2`. It adds two functions giving access
to a `double` marker, that is useful for the mesh optimizers to keep
the mesh density everywhere while modifying the mesh.

\cgalRefines `TriangulationVertexBase_2`

\cgalHasModel `CGAL::Delaunay_mesh_vertex_base_2<Traits, Vb>`


*/

class DelaunayMeshVertexBase_2 {
public:

/// \name Access Functions
/// @{

/*!
returns the size requested locally around the vertex
*/
const double& sizing_info() const;

/*!
sets the size requested locally around the vertex
*/
void set_sizing_info(const double& s);

/// @}

}; /* end DelaunayMeshVertexBase_2 */

