namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_cost_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the Garland-Heckbert collapse cost of top edge in the priority
queue is larger than a certain threshold.
This predicate is meant to be used with `GarlandHeckbert_cost`.

\tparam FT is the number type of the point coordinates.

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
  GarlandHeckbert_cost_stop_predicate(FT threshold);
  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns `current_cost >= threshold`.
  All other parameters are ignored (but exist since this is a generic policy).
  */
  template <typename F, typename Profile>
  bool operator()(const F& current_cost,
                  const Profile&,
                  std::size_t,
                  std::size_t) const;
  /// @}
};
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
