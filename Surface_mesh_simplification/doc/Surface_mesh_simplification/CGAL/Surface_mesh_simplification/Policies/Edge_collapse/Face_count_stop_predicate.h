namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

\cgalModels{StopPredicate}

The class `Face_count_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the number of current faces drops below a certain threshold.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\sa `CGAL::Surface_mesh_simplification::Edge_count_stop_predicate<TriangleMesh>`
\sa `CGAL::Surface_mesh_simplification::Face_count_ratio_stop_predicate<TriangleMesh>`
*/
template <typename TriangleMesh>
class Face_count_stop_predicate
{
public:

  /// \name Creation
  /// @{

  /*!
  Initializes the predicate establishing the `threshold` value.
  */
  Face_count_stop_predicate(const Edge_profile::edges_size_type threshold);

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns `true` if the current face count is strictly smaller than `threshold`, and `false` otherwise.
  All other parameters are ignored (but exist since this is a generic policy).
  */
  bool operator()(const Edge_profile::FT& current_cost,
                  const Edge_profile& edge_profile,
                  const Edge_profile::edges_size_type initial_edge_count,
                  const Edge_profile::edges_size_type current_edge_count) const;

  /// @}

};

} // namespace Surface_mesh_simplification
} // namespace CGAL
