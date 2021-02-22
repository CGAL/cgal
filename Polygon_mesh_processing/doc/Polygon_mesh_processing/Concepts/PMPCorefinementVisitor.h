/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPCorefinementVisitor` defines the requirements for the visitor
/// used in \link PMP_corefinement_grp corefinement-related functions \endlink to track
/// the creation of new faces and new edges.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Corefinement::Default_visitor`.

class PMPCorefinementVisitor{
public:
/// Mesh type
typedef unspecified_type Triangle_mesh;
/// Face descriptor type
typedef unspecified_type face_descriptor;
/// Halfedge descriptor type
typedef unspecified_type halfedge_descriptor;

/// @name Functions used by corefine() when faces are split
/// @{
  /// called before the triangulation of `f_split` in `tm`. Note that `f_split`
  /// will be one of the faces of the triangulation. Each subsequent call to
  /// `before_subface_created()`/`after_subface_created()` will correspond to
  /// the creation of a new face of the triangulation of `f_split`.
  void before_subface_creations(face_descriptor f_split, Triangle_mesh& tm);
  /// called when the triangulation of a face in `tm` is finished
  void after_subface_creations(Triangle_mesh& tm);
  /// called before creating a new triangle face in `tm` to triangulate the face passed to `before_subface_creations()`
  void before_subface_created(Triangle_mesh& tm);
  /// called after creating a new triangle face `f_new` in `tm` to triangulate the face passed to `before_subface_creations()`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_subface_created(face_descriptor f_new, Triangle_mesh& tm);
/// @}

/// @name Functions used by corefine() when edges are split
/// @{
  /// called before the edge of `h` in `tm` is split. Each subsequent call to
  /// `edge_split()` until the call to `after_edge_split()` will correspond to
  /// the split of that edge. If `edge_split(h_i, tm)` is called for `i=1` to `n`,
  /// `h_1`, `h_2`, ... ,`h_n`, `h` is the sequence of halfedges representing the
  /// edge split (with the same initial orientation)
  void before_edge_split(halfedge_descriptor h, TriangleMesh& tm);
  /// called when a new split is done. The target of `hnew` is a new split vertex.
  void edge_split(halfedge_descriptor hnew, TriangleMesh& tm);
  /// called when the split of the halfedge `h` passed at the later call to `before_edge_split()` is finished.
  void after_edge_split();
/// @}

/// @name Functions used by Boolean operations functions using corefinement.
/// These functions are not needed if you only call `corefine()`.
/// @{
  /// called before importing the face `f_src` of `tm_src` in `tm_tgt`
  void before_face_copy(face_descriptor f_src, Triangle_mesh& tm_src, Triangle_mesh& tm_tgt);
  /// called after importing the face `f_src` of `tm_src` in `tm_tgt`. The new face is `f_tgt`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_face_copy(face_descriptor  f_src, Triangle_mesh& tm_src,
                       face_descriptor  f_tgt, Triangle_mesh& tm_tgt);
/// @}
};
