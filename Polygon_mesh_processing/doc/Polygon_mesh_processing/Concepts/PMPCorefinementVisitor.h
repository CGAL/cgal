/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPCorefinementVisitor` defines the requirements for the visitor
/// used in \link PMP_corefinement_grp corefinement-related functions \endlink to track
/// the creation of new faces.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Corefinement::Default_visitor`.


class PMPCorefinementVisitor{
public:
/// Mesh type
typedef unspecified_type Triangle_mesh;
/// Face decriptor type
typedef unspecified_type face_descriptor;

/// @name Functions used by corefine()
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
