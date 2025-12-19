
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_site_2` is a model for the concept
`SegmentDelaunayGraphSite_2`.
\tparam K must be a model of the `Kernel` concept.

\cgalModels{SegmentDelaunayGraphSite_2}

\sa `Kernel`
\sa `SegmentDelaunayGraphSite_2`
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename K >
class Segment_Delaunay_graph_site_2 {
public:

/// \name Types
/// The class `Segment_Delaunay_graph_site_2` introduces the following
/// type in addition to the types in the concept
/// `SegmentDelaunayGraphSite_2`.
/// @{

/*!
A type for the template parameter `K`.
*/
typedef K Kernel;

/// @}

}; /* end Segment_Delaunay_graph_site_2 */
} /* end namespace CGAL */
