
namespace CGAL {

/*!
\ingroup PkgBasicViewerClasses

The class `Drawing_functor` is used to tune the way that a given data-structure of CGAL is drawn.
The different std::function can be modified to change the behavior of the drawing.

\tparam DS a data structure of CGAL.
\tparam VertexDescriptor a descriptor of vertices of DS.
\tparam EdgeDescriptor a descriptor of edges of DS.
\tparam FaceDescriptor a descriptor of faces of DS.

*/
template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor>
struct Drawing_functor
{
public:
  /// std::function that returns `true` if the given vertex must be drawn, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, VertexDescriptor)> draw_vertex;

  /// std::function that returns `true` if the given edge must be drawn, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, EdgeDescriptor)>   draw_edge;

  /// std::function that returns `true` if the given face must be drawn, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, FaceDescriptor)>   draw_face;

  /// std::function that returns `true` if the given vertex is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, VertexDescriptor)> colored_vertex;

  /// std::function that returns `true` if the given edge is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, EdgeDescriptor)>   colored_edge;

  /// std::function that returns `true` if the given face is colored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, FaceDescriptor)>   colored_face;

  /// std::function that returns `true` if the given face is drawn in wireframe, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, FaceDescriptor)> face_wireframe;

  /// std::function that returns the color of the given vertex.
  /// nullptr by default.
  std::function<CGAL::IO::Color(const DS &, VertexDescriptor)> vertex_color;

  /// std::function that returns the color of the given edge.
  /// nullptr by default.
  std::function<CGAL::IO::Color(const DS &, EdgeDescriptor)>   edge_color;

  /// std::function that returns the color of the given face.
  /// nullptr by default.
  std::function<CGAL::IO::Color(const DS &, FaceDescriptor)>   face_color;

  /// Disable the drawing of vertices.
  void disable_vertices();

  /// Enable the drawing of vertices.
  void enable_vertices();

  /// Disable the drawing of edges.
  void disable_edges();

  /// Enable the drawing of edges.
  void enable_edges();

    /// Disable the drawing of faces.
  void disable_faces();

  /// Enable the drawing of faces.
  void enable_faces();
};

/*!
\ingroup PkgBasicViewerClasses

The class `Drawing_functor_with_volume` is used to tune the way that a given data-structure of CGAL is drawn, for a data-structure that contains volumes.
The different std::function can be modified to change the behavior of the drawing.

\tparam DS a data structure of CGAL.
\tparam VertexDescriptor a descriptor of vertices of DS.
\tparam EdgeDescriptor a descriptor of edges of DS.
\tparam FaceDescriptor a descriptor of faces of DS.
\tparam VolumeDescriptor a descriptor of volumes of DS.

*/

// Drawing functor for a 3D data structure
// (with vertices, edges, faces and volumes)
template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor,
          typename VolumeDescriptor>
struct Drawing_functor_with_volume :
    public Drawing_functor<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>
{
public:
  /// std::function that returns `true` if the given volume must be drawn, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, VolumeDescriptor)> draw_volume;

  /// std::function that returns the color of the given volume.
  /// Returns `false` by default.
  std::function<bool(const DS &, VolumeDescriptor)> colored_volume;

  /// std::function that returns `true` if the given volume is drawn in wireframe, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, VolumeDescriptor)> volume_wireframe;

  /// std::function that returns the color of the given volume.
  /// nullptr by default.
  std::function<CGAL::IO::Color(const DS &, VolumeDescriptor)> volume_color;

  /// Disable the drawing of volumes.
  void disable_volumes();

  /// Enable the drawing of volumes.
  void enable_volumes();
};

} // End namespace CGAL
