
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Points

The class `Regular_triangulation_adaptation_traits_2` provides a model for the `AdaptationTraits_2`
concept. The template parameter of the `Regular_triangulation_adaptation_traits_2` class must be a
model of the `DelaunayGraph_2` concept, and in particular it has
the semantics of a 2D regular triangulation.

\cgalModels{AdaptationTraits_2}

\sa `AdaptationTraits_2`
\sa `DelaunayGraph_2`
\sa `Voronoi_diagram_2<DG,AT,AP>`
\sa `CGAL::Regular_triangulation_2<Traits,Tds>`

*/
template< typename RT2 >
struct Regular_triangulation_adaptation_traits_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Has_nearest_site_2;

/// @}

}; /* end Regular_triangulation_adaptation_traits_2 */
} /* end namespace CGAL */
