
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Edge_length_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the top edge in the priority queue is larger than a certain threshold.
This predicate is meant to be used with `Edge_length_cost`.

\tparam FT is the number type of the point coordinates, and must be equal to `Edge_profile::FT`.

\cgalModels{StopPredicate}

*/
template <typename FT>
class Edge_length_stop_predicate
{
public:

  /// \name Creation
  /// @{

  /*!
  Initializes the predicate establishing the `threshold` value.
  */
  Edge_length_stop_predicate<TriangleMesh>(const Edge_profile::FT threshold);

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns `(CGAL::squared_distance(edge_profile.p0(),edge_profile.p1()) > threshold*threshold)`.
  All other parameters are ignored (but exist since this is a generic policy).
  */
  bool operator()(const Edge_profile::FT&,
                  const Edge_profile& edge_profile,
                  const Edge_profile::edges_size_type initial_edge_count,
                  const Edge_profile::edges_size_type current_edge_count) const;

  /// @}

};

} // namespace Surface_mesh_simplification
} // namespace CGAL
