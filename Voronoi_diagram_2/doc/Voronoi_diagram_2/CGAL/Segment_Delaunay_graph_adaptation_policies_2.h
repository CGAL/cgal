
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Segments

The class `Segment_Delaunay_graph_caching_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Segment_Delaunay_graph_caching_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D segment Delaunay graph. This policy
caches the results of the edge and face rejectors and results in a
Voronoi diagram that has no degenerate features, i.e., no Voronoi
edges of zero length and no Voronoi faces of zero area.

\cgalModels `AdaptationPolicy_2`

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>`

*/
template< typename SDG2 >
struct Segment_Delaunay_graph_caching_degeneracy_removal_policy_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_false Has_inserter;

/// @}

}; /* end Segment_Delaunay_graph_caching_degeneracy_removal_policy_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Segments

The class `Segment_Delaunay_graph_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Segment_Delaunay_graph_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D segment Delaunay graphs. This policy
results in a Voronoi diagram that has no degenerate features,
i.e., it has no Voronoi edges of zero length and no Voronoi faces of
zero area.

\cgalModels `AdaptationPolicy_2`

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Segment_Delaunay_graph_caching_degeneracy_removal_policy_2<SDG2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>`

*/
template< typename SDG2 >
struct Segment_Delaunay_graph_degeneracy_removal_policy_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_inserter;

/// @}

}; /* end Segment_Delaunay_graph_degeneracy_removal_policy_2 */
} /* end namespace CGAL */
