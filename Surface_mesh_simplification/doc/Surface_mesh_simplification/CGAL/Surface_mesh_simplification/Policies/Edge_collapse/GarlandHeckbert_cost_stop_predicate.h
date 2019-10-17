namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_cost_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the Garland-Heckbert collapse cost of top edge in the priority
queue is larger than a certain threshold.
This predicate is meant to be used with `GarlandHeckbert_cost`.

\tparam FT is the number type of the point coordinates and must be equal to `Edge_profile::FT`.

\cgalModels `StopPredicate`

*/
template <class FT>
class GarlandHeckbert_cost_stop_predicate
{
public:
  /// \name Creation
  /// @{

  /*!
  Initializes the predicate establishing the `threshold` value.
  */
  GarlandHeckbert_cost_stop_predicate(Edge_profile::FT threshold);
  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns `current_cost >= threshold`.
  */
  bool operator()(const Edge_profile::FT& current_cost,
                  const Edge_profile&,
                  Edge_profile::edges_size_type initial_edge_count,
                  Edge_profile::edges_size_type current_edge_count) const;
  /// @}
};
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
