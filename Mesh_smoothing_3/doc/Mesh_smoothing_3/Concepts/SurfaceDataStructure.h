/*!
\ingroup PkgMeshSmoothing3Concepts
\cgalConcept

The concept `SurfaceDataStructure` describes the way the surface mesh will be accessed.

\sa `CGAL::Mesh_smoothing_3::Mesh_smoother`
\sa `MeshDataStructure`

*/
class SurfaceDataStructure {
public:

/// \name Types
/// @{

/*!
Descriptor used to access a face information
*/
using Face_descriptor = unspecified_type;


/*!
Vector type.
*/
using Normal_3 = unspecified_type;

/*!
Index associated with a surface patch to identify the patch it belongs to. This is used to query the patch information from the user.
*/
using Surface_patch_index = unspecified_type;

/// @}

/// \name Operations
/// The following functions are used to access surface data:
/// @{

/*!
Used only to reserve memory. std::size_t is optional but will avoid warnings.
*/
std::size_t nb_faces() const;

/*!
Provides an iterable range over the Face_descriptor of the mesh
*/
unspecified_type face_range() const;

/*!
Returns the number of vertices of given face. std::size_t is optional but will avoid warnings.
*/
std::size_t nb_face_vertices(Face_descriptor face) const;

/*!
Returns an identifier (patch id, face id, ...) related to the given face.
*/
Surface_patch_index patch_id(Face_descriptor face) const;

/*!
Provides an iterable range of Vertex_descriptor as defined in `MeshDataStructure` to iterate over the vertices of a face. 
*/
unspecified_type face_vertices(Face_descriptor face) const;

/// @}



}; /* end SurfaceDataStructure */
