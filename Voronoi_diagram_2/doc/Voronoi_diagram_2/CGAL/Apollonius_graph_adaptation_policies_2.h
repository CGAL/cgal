
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Disks

The class `Apollonius_graph_caching_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Apollonius_graph_caching_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D Apollonius graph. This policy
caches the results of the edge and face rejectors and results in a
Voronoi diagram that has no degenerate features, i.e., no Voronoi
edges of zero length and no Voronoi faces of zero area.

\cgalModels{AdaptationPolicy_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Apollonius_graph_degeneracy_removal_policy_2<AG2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>`

*/
template< typename AG2 >
struct Apollonius_graph_caching_degeneracy_removal_policy_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_inserter;

/// @}

}; /* end Apollonius_graph_caching_degeneracy_removal_policy_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Disks

The class `Apollonius_graph_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Apollonius_graph_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D Apollonius graph. This policy
results in a Voronoi diagram that has no degenerate features,
i.e., it has no Voronoi edges of zero length and no Voronoi faces of
zero area.

\cgalModels{AdaptationPolicy_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Apollonius_graph_caching_degeneracy_removal_policy_2<AG2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>`

*/
template< typename AG2 >
struct Apollonius_graph_degeneracy_removal_policy_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_inserter;

/// @}

}; /* end Apollonius_graph_degeneracy_removal_policy_2 */
} /* end namespace CGAL */
