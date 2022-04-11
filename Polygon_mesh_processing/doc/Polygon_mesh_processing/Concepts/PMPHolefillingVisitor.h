/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPHolefillingVisitor` defines the requirements for the visitor
/// used in \link PMP_hole_filling_grp hole-filling-related functions \endlink to track
/// the creation of new faces and new edges.
/// The hole filling may use a 2D constrained triangulation
/// for almost planar holes.  If that is not appropriate or fails it
/// may use an algorithm with at most quadratic running time. If that fails
/// it uses an algorithm with cubic running time.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Holefilling::Default_visitor`.

class PMPHolefillingVisitor{
public:

  /// called when the planar hole filling algorithm starts.
  void start_planar_phase() const;

  /// called when the  planar hole filling algorithm stops.
  /// @param success `true` when the hole could be filled in this phase.
  void end_planar_phase(bool success) const;

  /// called when the algorithm with quadratic running time starts.
  /// @param N the upperbound on the number of steps
  void start_quadratic_phase(int N) const;

    /// called at each step. There may be less than `N` calls as this
    /// is an upperbound.
  void quadratic_step() const;

  /// called when the algorithm with quadratic running time ends.
    /// @param success `true` when the hole could be filled in this
    /// phase.
  void end_quadratic_phase(bool success) const;

  /// called when the algorithm with cubic running time starts.
    /// @param N the upperbound on the number of steps
  void start_cubic_phase(int N) const;

  /// called at each step. This will be called `N` times.
  void cubic_step() const;

  /// called when the algorithm with cubic running time ends.
  void end_cubic_phase() const;


};
