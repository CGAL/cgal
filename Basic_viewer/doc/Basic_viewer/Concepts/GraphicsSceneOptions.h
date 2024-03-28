/*!
\ingroup PkgBasicViewerConcepts

The concept `GraphicsSceneOptions` defines data and methods used to tune the way that the cells of a given data structure of \cgal are considered for drawing or to be added into a graphics scene.
The different `std::function` can be modified to change for example the behavior of the drawing.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Graphics_scene_options `CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>`\endlink}
\cgalHasModelsEnd

*/
class GraphicsSceneOptions
{
public:
  /*!
    A data structure of \cgal.
  */
  typedef unspecified_type DS;

  /*!
    A descriptor of vertices of `DS`
  */
  typedef unspecified_type vertex_descriptor;

  /*!
    A descriptor of edges of `DS`
  */
  typedef unspecified_type edge_descriptor;

  /*!
    A descriptor of faces of `DS`
  */
 typedef unspecified_type face_descriptor;

  /// `std::function` that returns `true` if the given vertex must be ignored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, vertex_descriptor)> ignore_vertex;

  /// `std::function` that returns `true` if the given edge must be ignored, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, edge_descriptor)> ignore_edge;

  /// `std::function` that returns `true` if the given face must be ignored, `false` otherwise.
  /// Returns `true` by default.
  std::function<bool(const DS &, face_descriptor)> ignore_face;

  /// `std::function` that returns `true` if the given vertex is colored, `false` otherwise.
  /// Returns `false` by default.
  /// For non colored vertices, this is the role of the user of a graphic scene to decide which color must be used (cf. for example `Basic_viewer`, `vertices_mono_color`).
  std::function<bool(const DS &, vertex_descriptor)> is_vertex_colored;

  /// `std::function` that returns `true` if the given edge is colored, `false` otherwise.
  /// For non colored edges, this is the role of the user of a graphic scene to decide which color must be used (cf. for example `Basic_viewer`, `edges_mono_color`).
  /// Returns `false` by default.
  std::function<bool(const DS &, edge_descriptor)> is_edge_colored;

  /// `std::function` that returns `true` if the given face is colored, `false` otherwise.
  /// For non colored faces, this is the role of the user of a graphic scene to decide which color must be used (cf. for example `Basic_viewer`, `faces_mono_color`).
  /// Returns `false` by default.
  std::function<bool(const DS &, face_descriptor)> is_face_colored;

  /// `std::function` that returns `true` if the given face is in wireframe, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, face_descriptor)> is_face_wireframe;

  /// `std::function` that returns the color of the given vertex.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, vertex_descriptor)> vertex_color;

  /// `std::function` that returns the color of the given edge.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, edge_descriptor)> edge_color;

  /// `std::function` that returns the color of the given face.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, face_descriptor)> face_color;

  /// ignores all vertices when `b` is `true`; otherwise ignores only vertices for which `ignore_vertex()` returns `true`.
  void ignore_all_vertices(bool b);

  /// ignores all edges when `b` is `true`; otherwise ignores only edges for which `ignore_edge()` returns `true`.
  void ignore_all_edges(bool b);

  /// ignores all faces when `b` is `true`; otherwise ignores only faces for which `ignore_face()` returns `true`.
  void ignore_all_faces(bool b);
};

/*!
\ingroup PkgBasicViewerConcepts
\cgalConcept

The concept `GraphicsSceneOptionsWithVolumes` extends the concept `GraphicsSceneOptions` to deal with data structures that represent volumes.

\cgalRefines{GraphicsSceneOptions}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Graphics_scene_options `CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor, VolumeDescriptor>`\endlink}
\cgalHasModelsEnd

*/
class GraphicsSceneOptionsWithVolumes
{
public:
  /*!
    %A descriptor of volumes of `DS`.
  */
  typedef unspecified_type volume_descriptor;

  /// `std::function` that returns `true` if the given volume must be ignored, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, volume_descriptor)> ignore_volume;

  /// `std::function` that returns `true` if the given volume is colored, `false` otherwise.
  /// For non colored volumes, this is the role of the user of a graphic scene to decide which color must be used (cf. for example `Basic_viewer`, `faces_mono_color`).
  /// Returns `false` by default.
  std::function<bool(const DS &, volume_descriptor)> is_volume_colored;

  /// `std::function` that returns `true` if the given volume is in wireframe, `false` otherwise.
  /// Returns `false` by default.
  std::function<bool(const DS &, volume_descriptor)> is_volume_wireframe;

  /// `std::function` that returns the color of the given volume, i.e. the color of all the faces of this volume.
  /// `nullptr` by default.
  std::function<CGAL::IO::Color(const DS &, volume_descriptor)> volume_color;

  /// ignores all volumes when `b` is `true`; otherwise ignore only volumes for which `ignore_volume()` returns `true`.
  void ignore_all_volumes(bool b);
};
