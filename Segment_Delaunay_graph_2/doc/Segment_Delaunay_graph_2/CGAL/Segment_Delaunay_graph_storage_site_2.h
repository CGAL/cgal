
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_storage_site_2` is a model for the concept
`SegmentDelaunayGraphStorageSite_2`.

\tparam Gt must be a model of the `SegmentDelaunayGraphStorageTraits_2` concept.

\cgalModels{SegmentDelaunayGraphStorageSite_2}

\sa `CGAL::Segment_Delaunay_graph_site_2<K>`
*/
template< typename St >
class Segment_Delaunay_graph_storage_site_2 {
public:

/// \name Types
/// The class `Segment_Delaunay_graph_storage_site_2` introduces the
/// following type in addition to the types in the concept
/// `SegmentDelaunayGraphStorageSite_2`.
/// @{

/*!
A type for the template parameter `St`.
*/
typedef St Storage_traits;

/// @}

}; /* end Segment_Delaunay_graph_storage_site_2 */
} /* end namespace CGAL */
