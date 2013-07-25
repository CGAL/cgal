
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagramAdaptor2Points

The class `Delaunay_triangulation_adaptation_traits_2` provides a model for the `AdaptationTraits_2` 
concept. The template parameter of the `Delaunay_triangulation_adaptation_traits_2` class must be a 
model of the `DelaunayGraph_2` concept, and in particular it has 
the semantics of a 2D Delaunay triangulation. 

\cgalModels `AdaptationTraits_2`

\sa `AdaptationTraits_2` 
\sa `DelaunayGraph_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>` 

*/
template< typename DT2 >
class Delaunay_triangulation_adaptation_traits_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef CGAL::Tag_true Has_nearest_site_2; 

/// @}

}; /* end Delaunay_triangulation_adaptation_traits_2 */
} /* end namespace CGAL */
