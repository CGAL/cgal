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
  void before_subface_creations(face_descriptor f_split, const Triangle_mesh& tm);
  /// called when the triangulation of a face in `tm` is finished.
  void after_subface_creations(const Triangle_mesh& tm);
  /// called before creating a new triangle face in `tm` to triangulate the face passed to `before_subface_creations()`.
  void before_subface_created(const Triangle_mesh& tm);
  /// called after creating a new triangle face `f_new` in `tm` to triangulate the face passed to `before_subface_creations()`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_subface_created(face_descriptor f_new, const Triangle_mesh& tm);
/// @}

/// @name Functions used by corefine() when edges are split or created
/// @{
  /// called before the edge of `h` in `tm` is split. Each subsequent call to
  /// `edge_split()` until the call to `after_edge_split()` will correspond to
  /// the split of that edge. If `edge_split(h_i, tm)` is called for `i=1` to `n`,
  /// `h_1`, `h_2`, ... ,`h_n`, `h` is the sequence of halfedges representing the
  /// edge split (with the same initial orientation). There is only one call per edge.
  void before_edge_split(halfedge_descriptor h, const Triangle_mesh& tm);
  /// called when a new split is done. The target of `hnew` is a new split vertex. There is only one call per edge.
  void edge_split(halfedge_descriptor hnew, const Triangle_mesh& tm);
  /// called when the split of the halfedge `h` passed at the later call to `before_edge_split()` is finished.
  void after_edge_split();
  ///  called when a new edge has been added to triangulate a face. The face triangulated is `f_split`
  ///  in the last call to `before_subface_creations(f_split, tm)`. There is only one call per edge.
  void add_retriangulation_edge(halfedge_descriptor h, const Triangle_mesh& tm);
/// @}

/// @name Functions used by corefine() when a new vertex is created.
/// @{

  /// called when a new intersection point is detected. The intersection is detected using
  /// a face of `tm_f` and an edge of `tm_e`.
  /// \param i_id the id of the intersection point, starting at 0. Ids are consecutive.
  /// \param sdim indicates the dimension of the simplex part of the face that is intersected by the edge
  ///             (0 for a vertex, 1 for an edge, 2 for the interior of the face)
  /// \param h_f a halfedge from `tm_f` indicating the simplex intersected:
  ///            if `sdim==0` the target of `h_f` is the intersection point,
  ///            if `sdim==1` the edge of `h_f` contains the intersection point in its interior,
  ///            if `sdim==2` the face of `h_f` contains the intersection point in its interior.
  /// \param h_e a halfedge from `tm_e`
  /// \param is_target_coplanar `true` iff the target of `h_e` is the intersection point
  /// \param is_source_coplanar `true` iff the source of `h_e` is the intersection point
  /// \param tm_f mesh containing `h_f`
  /// \param tm_e mesh containing `h_e`
  void intersection_point_detected(std::size_t i_id,
                                   int sdim,
                                   halfedge_descriptor h_f,
                                   halfedge_descriptor h_e,
                                   const Triangle_mesh& tm_f,
                                   const Triangle_mesh& tm_e,
                                   bool is_target_coplanar,
                                   bool is_source_coplanar);

  /// called when a new vertex is added in `tm` (either an edge split or a vertex inserted in the interior of a face).
  /// `i_id` is the intersection point id reported in `new_node_added`.
  /// For each mesh, a vertex with a given id will be reported exactly once, except if it is already an existing vertex.
  void new_vertex_added(std::size_t i_id, vertex_descriptor v, const Triangle_mesh& tm);
/// @}

/// @name Functions used by Boolean operations functions using corefinement.
/// These functions are not needed if you only call `corefine()`.
/// @{
  /// called before importing the face `f_src` of `tm_src` in `tm_tgt`
  void before_face_copy(face_descriptor f_src, const Triangle_mesh& tm_src, const Triangle_mesh& tm_tgt);
  /// called after importing the face `f_src` of `tm_src` in `tm_tgt`. The new face is `f_tgt`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_face_copy(face_descriptor  f_src, const Triangle_mesh& tm_src,
                       face_descriptor  f_tgt, const Triangle_mesh& tm_tgt);
  /// called before importing the edge of `h_src` of `tm_src` in `tm_tgt`. There is one call per edge.
  void before_edge_copy(halfedge_descriptor h_src/, const Triangle_mesh& tm_src, const Triangle_mesh& tm_tgt){}
  /// called after importing the edge of `h_src` of `tm_src` in `tm_tgt`. The corresponding new halfedge is `h_tgt`.
  /// There is only one call per edge.
  void after_edge_copy(halfedge_descriptor h_src, const Triangle_mesh& tm_src,
                       halfedge_descriptor h_tgt, const Triangle_mesh& tm_tgt);
  /// called before a patch boundary edge is duplicated to disconnect patches of `tm`
  /// (When an in-place operation and an out-of-place are both requested).
  void before_edge_duplicated(halfedge_descriptor h, const Triangle_mesh& tm);
  /// called when the edge of `h_src` has been duplicated into `h_new` in `tm`
  /// (When an in-place operation and an out-of-place are both requested).
  void after_edge_duplicated(halfedge_descriptor h_src,
                             halfedge_descriptor h_new, const Triangle_mesh& tm);
  /// called when an intersection edge (represented in input meshes `tm_src1` and `tm_src2` by `h_src1` and `h_src2`,
  /// respectively) is imported in `tm_tgt` as `h_tgt`. There is only one call per edge.
  /// (Called only when an out-of-place operation is requested)
  void intersection_edge_copy(halfedge_descriptor h_src1, const Triangle_mesh& tm_src1,
                              halfedge_descriptor h_src2, const Triangle_mesh& tm_src2,
                              halfedge_descriptor h_tgt,  const Triangle_mesh& tm_tgt);
  /// called before vertex `v_src` from `tm_src` is copied in `tm_tgt`
  void before_vertex_copy(vertex_descriptor v_src, const Triangle_mesh& tm_src, const Triangle_mesh& tm_tgt);
  /// called after vertex `v_src` from `tm_src` is copied in `tm_tgt`. The new vertex is `v_tgt`.
  /// The point has already been put in the vertex point map.
  void after_vertex_copy(vertex_descriptor v_src, const Triangle_mesh& tm_src,
                         vertex_descriptor v_tgt, const Triangle_mesh& tm_tgt);
/// @}
};
