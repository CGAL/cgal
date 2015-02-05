
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2

The class `Segment_Delaunay_graph_storage_site_2` is a model for the concept 
`SegmentDelaunayGraphStorageSite_2`. 

\tparam Gt must be a model of the `SegmentDelaunayGraphTraits_2` concept. 

\cgalModels `SegmentDelaunayGraphStorageSite_2`

\sa `SegmentDelaunayGraphSite_2` 
\sa `SegmentDelaunayGraphTraits_2` 
\sa `CGAL::Segment_Delaunay_graph_site_2<K>` 
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>` 

*/
template< typename Gt >
class Segment_Delaunay_graph_storage_site_2 {
public:

/// \name Types 
/// The class `Segment_Delaunay_graph_storage_site_2` introduces the
/// following type in addition to the types in the concept
/// `SegmentDelaunayGraphStorageSite_2`.
/// @{

/*!
A type for the template parameter `Gt`. 
*/ 
typedef Gt Geom_traits; 

/// @}

}; /* end Segment_Delaunay_graph_storage_site_2 */
} /* end namespace CGAL */
