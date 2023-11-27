/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `EdgeCollapseSimplificationVisitor` describes the requirements for the <I>visitor object</I> which is used to track the edge collapse simplification algorithm.

The several callbacks given as member functions in the visitor are called from certain places in the algorithm implementation.

*/
class EdgeCollapseSimplificationVisitor {
public:

/// \name Types
/// @{

/// Convenience typedef to the Edge_profile class.
typedef CGAL::Surface_mesh_simplification::Edge_profile Edge_profile;

/*!
The type of the surface mesh to simplify. Must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
*/
typedef Edge_profile::Triangle_mesh TriangleMesh;

/// \name Operations
/// @{

/*!
Called before the algorithm starts.
*/
void OnStarted(TriangleMesh& surface_mesh);

/*!
Called after the algorithm finishes.
*/
void OnFinished(TriangleMesh& surface_mesh);

/*!
Called when the `StopPredicate` returned `true`
(but not if the algorithm terminates because the surface mesh could not be simplified any further).

*/
void OnStopConditionReached(TriangleMesh& surface_mesh);

/*!
Called during the <I>collecting phase</I> (when a cost is assigned to the edges),
for each edge collected.

*/
void OnCollected(const Edge_profile& profile,
                 std::optional<Edge_profile::FT> cost);

/*!
Called during the <I>processing phase</I> (when edges are collapsed),
for each edge that is selected.

This method is called before the algorithm checks if the edge is collapsible.

`cost` indicates the current collapse cost for the edge.
If absent (meaning that it could not be computed)
the edge will not be collapsed.

`initial_edge_count` and `current_edge_count` refer to the number of edges.

*/
void OnSelected(const Edge_profile& profile,
                std::optional<Edge_profile::FT> cost,
                const Edge_profile::edges_size_type initial_edge_count,
                const Edge_profile::edges_size_type current_edge_count);

/*!
Called when an edge is about to be collapsed and replaced by a vertex
whose position is `*placement`.

If `placement` is absent (meaning that it could not be computed)
the edge will not be collapsed.

*/
void OnCollapsing(const Edge_profile& profile,
                  std::optional<Edge_profile::Point> placement);

/*!
Called when an edge has been collapsed and replaced by the vertex `vd`
*/
void OnCollapsed(const Edge_profile& profile,
                 const Edge_profile::vertex_descriptor vd) {}

/*!
Called for each selected edge which cannot be
collapsed because doing so would change the topological
type of the surface mesh (turn it into a non-manifold
for instance).

*/
void OnNonCollapsable(const Edge_profile& profile);

/// @}

}; /* end EdgeCollapseSimplificationVisitor */
