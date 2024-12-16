
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

\cgalModels{StopPredicate}

The class `Edge_count_ratio_stop_predicate` is a model for the `StopPredicate` concept
which returns `true` when the relation between the initial and current number of edges drops below a certain ratio.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\sa `CGAL::Surface_mesh_simplification::Edge_count_stop_predicate<TriangleMesh>`
\sa `CGAL::Surface_mesh_simplification::Face_count_ratio_stop_predicate<TriangleMesh>`
*/
template< typename TriangleMesh>
class Edge_count_ratio_stop_predicate
{
public:

  /// \name Creation
  /// @{

  /*!
  Initializes the predicate establishing the `ratio`.
  */
  Edge_count_ratio_stop_predicate(const double ratio);

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns `(double(current_edge_count) / double(initial_edge_count)) < ratio`.
  All other parameters are ignored (but exist since this is a generic policy).
  */
  bool operator()(const Edge_profile::FT current_cost,
                  const Edge_profile& edge_profile,
                  const Edge_profile::edges_size_type initial_edge_count,
                  const Edge_profile::edges_size_type current_edge_count) const;

  /// @}

};

} // namespace Surface_mesh_simplification
} // namespace CGAL
