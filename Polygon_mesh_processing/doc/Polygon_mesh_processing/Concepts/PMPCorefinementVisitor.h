/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPCorefinementVisitor` defines the requirements for the visitor
/// used in \link PMP_corefinement_grp corefinement-related functions \endlink to track
/// the creation of new faces and new edges.
///
/// \cgalRefines{CopyConstructible}
/// \cgalHasModelsBegin
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Corefinement::Default_visitor}
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Corefinement::Non_manifold_output_visitor}
/// \cgalHasModelsEnd

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
  /// \param h_e a halfedge from `tm_e`
  /// \param h_f a halfedge from `tm_f` indicating the simplex intersected:
  ///            if `sdim==0` the target of `h_f` is the intersection point,
  ///            if `sdim==1` the edge of `h_f` contains the intersection point in its interior,
  ///            if `sdim==2` the face of `h_f` contains the intersection point in its interior.
  /// \param tm_e mesh containing `h_e`
  /// \param tm_f mesh containing `h_f`
  /// \param is_target_coplanar `true` iff the target of `h_e` is the intersection point
  /// \param is_source_coplanar `true` iff the source of `h_e` is the intersection point
  void intersection_point_detected(std::size_t i_id,
                                   int sdim,
                                   halfedge_descriptor h_e,
                                   halfedge_descriptor h_f,
                                   const Triangle_mesh& tm_e,
                                   const Triangle_mesh& tm_f,
                                   bool is_target_coplanar,
                                   bool is_source_coplanar);

  /// called when a new vertex is added in `tm` (either an edge split or a vertex inserted in the interior of a face).
  /// `i_id` is the intersection point id reported in `new_node_added`.
  /// For each mesh, a vertex with a given id will be reported exactly once, except if it is already an existing vertex.
  void new_vertex_added(std::size_t i_id, vertex_descriptor v, const Triangle_mesh& tm);
/// @}

/// @name Functions used by Boolean operations functions using corefinement.
/// These functions are not needed if only `corefine()` is called.
/// @{
  /// called before importing the face `f_src` of `tm_src` in `tm_tgt`.
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
  /// (Called only when an out-of-place operation is requested).
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

/// @name Functions used by corefine() for progress tracking
/// @{
  /// called before starting to detect intersections between faces of one mesh and edges of the other.
  void start_filtering_intersections();
  /// called during detection of intersections between faces of one mesh and edges of the other.
  /// `d` is a double value in `[0,1]` that is increasing with the number of calls. The closer
  /// `d`is to `1`, the closer the intersection detection is to completion.
  void progress_filtering_intersections(double d);
  /// called after detection of intersections between faces of one mesh and edges of the other.
  void end_filtering_intersections();

  /// called before processing intersections between the `n` pairs of coplanar faces.
  void start_handling_intersection_of_coplanar_faces(std::size_t n);
  /// called each time a pair of coplanar faces is processed.
  void intersection_of_coplanar_faces_step() const;
  /// called after processing all intersections between coplanar faces.
  void end_handling_intersection_of_coplanar_faces() const;

  /// called before processing intersections between edges and faces of two meshes (called twice).
  /// `n` is the number of edges possibly intersecting faces that will be processed.
  void start_handling_edge_face_intersections(std::size_t n);
  /// called each time an edge is processed.
  void edge_face_intersections_step();
  /// called after having processed edge-face intersections between two meshes.
  void end_handling_edge_face_intersections();

  /// called before triangulating the `n` split faces.
  void start_triangulating_faces(std::size_t n);
  /// called when triangulating one split face.
  void triangulating_faces_step();
  /// called after the triangulation of the split faces.
  void end_triangulating_faces();
/// @}

/// @name Functions used by Boolean operations functions using corefinement for progress tracking.
/// These functions are not needed if only `corefine()` is called.
/// called before computing the output of the Boolean operations after corefinement.
/// @{
  void start_building_output();
  /// called when the output of the Boolean operations is computed.
  void end_building_output();
  /// called before filtering intersection edges in the interior of a set of coplanar faces.
  void filter_coplanar_edges();
  /// called before segmenting input meshes in patches defined by connected components separated by intersection edges.
  void detect_patches();
  /// called before classifying which patches contribute to each Boolean operation.
  void classify_patches();
  /// called before classifying patches of `tm` that are free from intersection with the other mesh.
  void classify_intersection_free_patches(const TriangleMesh& tm);
  /// called before creating a new mesh for a Boolean operation of type `t`.
  void out_of_place_operation(Boolean_operation_type t);
  /// called before updating an input mesh to store the Boolean operation of type `t`.
  void in_place_operation(Boolean_operation_type t);
  /// called before updating both input meshes to store the Boolean operations of type `t1` and `t2`.
  void in_place_operations(Boolean_operation_type t1,Boolean_operation_type t2);
/// @}
};
