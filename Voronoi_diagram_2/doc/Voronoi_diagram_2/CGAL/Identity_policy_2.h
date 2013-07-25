
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagramAdaptor2

The class `Identity_policy_2` provides a model for the `AdaptationPolicy_2` 
concept. The first template parameter of the `Identity_policy_2` class must be a 
model of the `DelaunayGraph_2` concept, whereas as the second 
template parameter must be a model of the `AdaptationTraits_2` concept. 
This policy rejects no edge and no face of the Delaunay graph, thus 
giving a Voronoi diagram which is the true dual of the triangulation 
Delaunay graph. The Voronoi diagram created with this adaptation 
policy may have degenerate features, such as Voronoi edges of zero 
length, or Voronoi faces of zero area. This policy assumes that the 
Delaunay graph, that is adapted, allows for site insertions through an 
`insert` method that takes as argument an object of type 
`AT::Site_2`. The site inserter functor provided by this policy 
uses the aforementioned `insert` method. 

\cgalModels `AdaptationPolicy_2`

\sa `AdaptationTraits_2` 
\sa `DelaunayGraph_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 

*/
template< typename DG, typename AT >
class Identity_policy_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef CGAL::Tag_true Has_inserter; 

/// @}

}; /* end Identity_policy_2 */
} /* end namespace CGAL */
