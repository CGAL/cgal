
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagramAdaptor2Disks

The class `Apollonius_graph_adaptation_traits_2` provides a model for the `AdaptationTraits_2` 
concept. The template parameter of the `Apollonius_graph_adaptation_traits_2` class must be a 
model of the `DelaunayGraph_2` concept, and in particular it has 
the semantics of a (triangulated) 2D Apollonius graph. 

\cgalModels `AdaptationTraits_2`

\sa `AdaptationTraits_2` 
\sa `DelaunayGraph_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 
\sa `CGAL::Apollonius_graph_2<Gt,Agds>` 
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>` 

*/
template< typename AG2 >
class Apollonius_graph_adaptation_traits_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef CGAL::Tag_true Has_nearest_site_2; 

/// @}

}; /* end Apollonius_graph_adaptation_traits_2 */
} /* end namespace CGAL */
