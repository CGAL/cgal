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

/// \name Types
/// @{

/*!
The type of the surface mesh to simplify. Must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
*/
typedef unspecified_type TriangleMesh;

/*!
A field type representing the collapse cost
*/
typedef unspecified_type FT;

/*!
An integer type representing the number of edges
*/
typedef unspecified_type size_type;

/*!
The type of the edge profile cache. Must be a model of the `EdgeProfile` concept.
*/
typedef unspecified_type Profile;

/// @}

/// \name Operations
/// @{

/*!

This predicate is called each time an edge is selected for processing,
before it is collapsed.

`current_cost` is the cost of the selected edge.

`initial_count` and `current_count` are the number of initial and current edges.

If the return value is `true` the simplification terminates before processing the edge,
otherwise it continues normally.

*/
  bool operator()(const FT& current_cost,
                  const Profile& profile,
                  size_type initial_count,
                  size_type current_count) const;

/// @}

}; /* end StopPredicate */
