/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `StopPredicate` describes the requirements for the predicate which indicates if the simplification process must finish.

\cgalHasModel `CGAL::Surface_mesh_simplification::Count_stop_predicate<TriangleMesh>`
\cgalHasModel `CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<TriangleMesh>`
\cgalHasModel `CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<FT>`

*/
class StopPredicate {
public:

/// \name Operations
/// @{

/*!

This predicate is called each time an edge is selected for processing, before it is collapsed.
`current_edge_cost` is the cost of the selected edge.
`initial_edge_count` and `current_edge_count` are the number of initial and current edges.
If the return value is `true` the simplification terminates before processing the edge,
otherwise it continues normally.

*/
  bool operator()(const Edge_profile::FT& current_cost,
                  const Edge_profile& profile,
                  const Edge_profile::edges_size_type initial_edge_count,
                  const Edge_profile::edges_size_type current_edge_count) const;

/// @}

}; /* end StopPredicate */
