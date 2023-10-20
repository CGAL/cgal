
namespace CGAL {

/*!
\ingroup PkgBasicViewerClasses

The class `Graphics_scene_options` is used to tune the way that the cells of a given data structure of \cgal are considered.
The different `std::function` can be modified to change for example the behavior of the drawing.
`VolumeDescriptor` can be `void` for data structures that do not represent volumes. In such a case, all methods about volumes do not exist.

\tparam DS a data structure of \cgal.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.
\tparam VolumeDescriptor a descriptor of volumes of `DS`. `void` by default.

*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor,
          typename VolumeDescriptor=void>
struct Graphics_scene_options
{
public:
  /// `std::function` that returns `true` if the given vertex must be ignored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, VertexDescriptor)> ignore_vertex;

  /// `std::function` that returns `true` if the given edge must be ignored, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, EdgeDescriptor)> ignore_edge;

  /// `std::function` that returns `true` if the given face must be ignored, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, FaceDescriptor)> ignore_face;

  /// `std::function` that returns `true` if the given volume must be ignored, `false` otherwise.
  /// Exists only if `VolumeDescriptor` is not `void`.
  /// Returns `false` by default.
  std::function<bool(const DS &, VolumeDescriptor)> ignore_volume;

  /// `std::function` that returns `true` if the given vertex is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, VertexDescriptor)> is_vertex_colored;

  /// `std::function` that returns `true` if the given edge is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, EdgeDescriptor)> is_edge_colored;

  /// `std::function` that returns `true` if the given face is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, FaceDescriptor)> is_face_colored;

  /// `std::function` that returns `true` if the given volume is colored, `false` otherwise.
  /// Returns `false` by default.
  /// Exists only if `VolumeDescriptor` is not `void`.
  std::function<bool(const DS &, VolumeDescriptor)> is_volume_colored;

  /// `std::function` that returns `true` if the given face is in wireframe, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, FaceDescriptor)> is_face_wireframe;

  /// `std::function` that returns `true` if the given volume is in wireframe, `false` otherwise.
  /// Returns `false` by default.
  /// Exists only if `VolumeDescriptor` is not `void`.
  std::function<bool(const DS &, VolumeDescriptor)> is_volume_wireframe;

  /// `std::function` that returns the color of the given vertex.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, VertexDescriptor)> vertex_color;

  /// `std::function` that returns the color of the given edge.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, EdgeDescriptor)> edge_color;

  /// `std::function` that returns the color of the given face.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, FaceDescriptor)> face_color;

  /// `std::function` that returns the color of the given volume.
  /// `nullptr` by default.
  /// Exists only if `VolumeDescriptor` is not `void`.
  std::function<CGAL::IO::Color(const DS &, VolumeDescriptor)> volume_color;

  /// ignores all vertices when `b` is `true`; otherwise ignores only vertices for which `ignore_vertex()` returns `true`.
  void ignore_all_vertices(bool b);

  /// ignores all edges when `b` is `true`; otherwise ignores only edges for which `ignore_edge()` returns `true`.
  void ignore_all_edges(bool b);

  /// ignores all faces when `b` is `true`; otherwise ignores only faces for which `ignore_face()` returns `true`.
  void ignore_all_faces(bool b);

  /// ignores all volumes when `b` is `true`; otherwise ignore only volumes for which `ignore_volume()` returns `true`.
  /// Exists only if `VolumeDescriptor` is not `void`.
  void ignore_all_volumes(bool b);
};

} // End namespace CGAL
