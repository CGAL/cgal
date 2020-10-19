/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPTriangulateFaceVisitor` defines the requirements for the visitor
/// used in \link PMP_meshing_grp triangulation-related functions \endlink to track
/// the creation of new faces.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor`.


class PMPTriangulateFaceVisitor {
public:
/// Face decriptor type
typedef unspecified_type face_descriptor;

/// @name Functions used by triangulate_face() and triangulate_faces()
/// @{
  /// called before the triangulation of `f_split`. Note that `f_split`
  /// will be one of the faces of the triangulation. Each subsequent call to
  /// `after_subface_created()` will correspond to
  /// the creation of a new face of the triangulation of `f_split`.
  void before_subface_creations(face_descriptor f_split);

  /// called when the triangulation of a face in `tm` is finished.
  void after_subface_creations();

  /// called after creating a new triangle face `f_new` to triangulate the face passed to `before_subface_creations()`.
  /// Note that the call is placed just after a call to `add_face()` so the halfedge pointer is not set yet.
  void after_subface_created(face_descriptor f_new);

  /// @}

};
