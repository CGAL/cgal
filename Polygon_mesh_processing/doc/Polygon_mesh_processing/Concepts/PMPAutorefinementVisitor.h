/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPAutorefinementVisitor` defines the requirements for the visitor
/// used in `CGAL::Polygon_mesh_processing::autorefine_triangle_soup()` to track
/// the creation of new triangles.
///
/// \cgalRefines{CopyConstructible}
/// \cgalHasModelsBegin
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor}
/// \cgalHasModelsEnd

class PMPAutorefinementVisitor{
public:

/// @name Functions called only if at least one intersection has been found
/// @{
  /// called when the final number of output triangles is known, `nbt` being the total number of triangles in the output.
  void number_of_output_triangles(std::size_t nbt);
  /// called for triangle with no intersection, `tgt_id` is the position in the triangle container after calling
  /// `autorefine_triangle_soup()`, while `src_id` was its position before calling the function.
  void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id);
  /// called for each subtriangle created from a triangle with intersection, `tgt_id` is the position in the triangle container after calling
  /// `autorefine_triangle_soup()` of the subtriangle, while `src_id` was the position of the original support triangle before calling the function.
  void new_subtriangle(std::size_t tgt_id, std::size_t src_id);
/// @}
};
