
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Segments

The class `Segment_Delaunay_graph_adaptation_traits_2` provides a model for the `AdaptationTraits_2`
concept. The template parameter of the `Segment_Delaunay_graph_adaptation_traits_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of the 2D (triangulated) segment Delaunay graph.

\cgalModels{AdaptationTraits_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>`

*/
template< typename SDG2 >
struct Segment_Delaunay_graph_adaptation_traits_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_nearest_site_2;

/// @}

}; /* end Segment_Delaunay_graph_adaptation_traits_2 */
} /* end namespace CGAL */
