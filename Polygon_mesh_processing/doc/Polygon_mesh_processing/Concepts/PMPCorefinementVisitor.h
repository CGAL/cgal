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
  /// called when the triangulation of a face in `tm` is finished.
  void after_subface_creations(Triangle_mesh& tm);
  /// called before creating a new triangle face in `tm` to triangulate the face passed to `before_subface_creations()`.
  void before_subface_created(Triangle_mesh& tm);
  /// called after creating a new triangle face `f_new` in `tm` to triangulate the face passed to `before_subface_creations()`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_subface_created(face_descriptor f_new, Triangle_mesh& tm);
/// @}

/// @name Functions used by corefine() when edges are split or created
/// @{
  /// called before the edge of `h` in `tm` is split. Each subsequent call to
  /// `edge_split()` until the call to `after_edge_split()` will correspond to
  /// the split of that edge. If `edge_split(h_i, tm)` is called for `i=1` to `n`,
  /// `h_1`, `h_2`, ... ,`h_n`, `h` is the sequence of halfedges representing the
  /// edge split (with the same initial orientation). There is only one call per edge.
  void before_edge_split(halfedge_descriptor h, TriangleMesh& tm);
  /// called when a new split is done. The target of `hnew` is a new split vertex. There is only one call per edge.
  void edge_split(halfedge_descriptor hnew, TriangleMesh& tm);
  /// called when the split of the halfedge `h` passed at the later call to `before_edge_split()` is finished.
  void after_edge_split();
  ///  called when a new edge has been added to triangulate a face. The face triangulated is `f_split`
  ///  in the last call to `before_subface_creations(f_split, tm)`. There is only one call per edge.
  void add_retriangulation_edge(halfedge_descriptor h, TriangleMesh& tm);
/// @}

/// @name Functions used by Boolean operations functions using corefinement.
/// These functions are not needed if you only call `corefine()`.
/// @{
  /// called before importing the face `f_src` of `tm_src` in `tm_tgt`
  void before_face_copy(face_descriptor f_src, const Triangle_mesh& tm_src, Triangle_mesh& tm_tgt);
  /// called after importing the face `f_src` of `tm_src` in `tm_tgt`. The new face is `f_tgt`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_face_copy(face_descriptor  f_src, const Triangle_mesh& tm_src,
                       face_descriptor  f_tgt, Triangle_mesh& tm_tgt);
  /// called before importing the edge of `h_src` of `tm_src` in `tm_tgt`. There is one call per edge.
  void before_edge_copy(halfedge_descriptor h_src/, const TriangleMesh& tm_src, TriangleMesh& tm_tgt){}
  /// called after importing the edge of `h_src` of `tm_src` in `tm_tgt`. The corresponding new halfedge is `h_tgt`.
  /// There is only one call per edge.
  void after_edge_copy(halfedge_descriptor h_src, const TriangleMesh& tm_src,
                       halfedge_descriptor h_tgt, TriangleMesh& tm_tgt);
  /// called before a patch boundary edge is duplicated to disconnect patches of `tm`
  /// (When an in-place operation and an out-of-place are both requested).
  void before_edge_duplicated(halfedge_descriptor h, TriangleMesh& tm);
  /// called when the edge of `h_src` has been duplicated into `h_new` in `tm`
  /// (When an in-place operation and an out-of-place are both requested).
  void after_edge_duplicated(halfedge_descriptor h_src,
                             halfedge_descriptor h_new, TriangleMesh& tm);
  /// called when an intersection edge (represented in input meshes `tm_src1` and `tm_src2` by `h_src1` and `h_src2`,
  /// respectively) is imported in `tm_tgt` as `h_tgt`. There is only one call per edge.
  /// (Called only when an out-of-place operation is requested)
  void intersection_edge_copy(halfedge_descriptor h_src1, const TriangleMesh& tm_src1,
                              halfedge_descriptor h_src2, const TriangleMesh& tm_src2,
                              halfedge_descriptor h_tgt,  TriangleMesh& tm_tgt);
/// @}
};
