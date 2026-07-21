/*!
\ingroup PkgMeshSmoothing3Concepts
\cgalConcept

The concept `PolylinesDataStructure` describes the way the curves on the mesh will be accessed.

\sa `CGAL::Mesh_smoothing_3::Mesh_smoother`
\sa `MeshDataStructure`
\sa `SurfaceDataStructure`

*/
class PolylinesDataStructure {
public:

/// \name Types
/// @{

/*!
Descriptor used to access a edge information
*/
using Edge_descriptor = unspecified_type;

/*!
Index associated with an edge to identify the curve it belongs to. This is used to query the curve information from the user.
*/
using Curve_index = unspecified_type;

/// @}

/// \name Operations
/// The following functions are used to access surface data:
/// @{

/*!
Used only to reserve memory. std::size_t is optional but will avoid warnings.
*/
std::size_t nb_edges() const;

/*!
Provides an iterable range over the Edge_descriptor of the mesh
*/
unspecified_type edge_range() const;

/*!
Returns an identifier (curve id, segment id, ...) related to the given edge.
*/
Curve_index curve_id(Edge_descriptor edge) const;

/*!
Return the ith vertex of the edge (max 2). Vertex_descriptor as given in `MeshDataStructure`. 
*/
Vertex_descriptor edge_vertex(Edge_descriptor edge, unsigned i) const;

/// @}


}; /* end PolylinesDataStructure */
