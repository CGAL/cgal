/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPHolefillingVisitor` defines the requirements for the visitor
/// used in \link PMP_hole_filling_grp hole-filling-related functions \endlink.
/// The hole filling may use a 2D constrained triangulation
/// for almost planar holes (*planar phase*).  If that is not appropriate or fails it
/// may use an algorithm with a quadratic running time relying on the 3D Delaunay triangulation (*quadratic phase*).
/// If that fails, it uses an algorithm with cubic running time (*cubic phase*).
///
/// \cgalRefines{CopyConstructible}
/// \cgalHasModelsBegin
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Hole_filling::Default_visitor}
/// \cgalHasModelsEnd

class PMPHolefillingVisitor{
public:

  /// called when the planar hole filling algorithm starts.
  void start_planar_phase() const;

  /// called when the  planar hole filling algorithm stops.
  /// @param success `true` when the hole could be filled in this phase.
  void end_planar_phase(bool success) const;

  /// called when the algorithm with quadratic running time starts.
  /// @param n the upperbound on the number of steps
  void start_quadratic_phase(std::size_t n) const;

    /// called at each step. There may be less than `n` calls as this
    /// is an upperbound.
  void quadratic_step() const;

  /// called when the algorithm with quadratic running time ends.
    /// @param success `true` when the hole could be filled in this
    /// phase.
  void end_quadratic_phase(bool success) const;

  /// called when the algorithm with cubic running time starts.
    /// @param n the upperbound on the number of steps
  void start_cubic_phase(std::size_t n) const;

  /// called at each step. This will be called `n` times.
  void cubic_step() const;

  /// called when the algorithm with cubic running time ends.
  void end_cubic_phase() const;

  /// called before refining the triangulation of the hole
  ///(`triangulate_and_refine_hole()` and `triangulate_refine_and_fair_hole()` only)
  void start_refine_phase();

  /// called after having refined the triangulation of the hole
  ///(`triangulate_and_refine_hole()` and `triangulate_refine_and_fair_hole()` only)
  void end_refine_phase();

  /// called before fairing the triangulation of the hole
  ///(`triangulate_refine_and_fair_hole()` only)
  void start_fair_phase();

  /// called after having faired the triangulation of the hole
  ///(`triangulate_refine_and_fair_hole()` only)
  void end_fair_phase();

};
