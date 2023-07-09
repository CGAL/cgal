
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Points

The class `Regular_triangulation_caching_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Regular_triangulation_caching_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D regular triangulation. This policy
caches the results of the edge and face rejectors and results in a
Voronoi diagram that has no degenerate features, i.e., no Voronoi
edges of zero length.

\cgalModels{AdaptationPolicy_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Regular_triangulation_2<Traits,Tds>`

*/
template< typename RT2 >
struct Regular_triangulation_caching_degeneracy_removal_policy_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_inserter;

/// @}

}; /* end Regular_triangulation_caching_degeneracy_removal_policy_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Points

The class `Regular_triangulation_degeneracy_removal_policy_2` provides a model for the `AdaptationPolicy_2`
concept. The template parameter of the `Regular_triangulation_degeneracy_removal_policy_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a (triangulated) 2D regular triangulation. This policy
results in a power diagram that has no degenerate features,
i.e., it has no Voronoi edges of zero length.

\cgalModels{AdaptationPolicy_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Regular_triangulation_caching_degeneracy_removal_policy_2<RT2>`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Regular_triangulation_2<Traits,Tds>`

*/
template< typename RT2 >
struct Regular_triangulation_degeneracy_removal_policy_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_inserter;

/// @}

}; /* end Regular_triangulation_degeneracy_removal_policy_2 */
} /* end namespace CGAL */
