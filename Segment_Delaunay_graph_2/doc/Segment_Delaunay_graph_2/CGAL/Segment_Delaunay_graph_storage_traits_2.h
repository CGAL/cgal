namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_storage_traits_2` provides a model for the
`SegmentDelaunayGraphStorageTraits_2` concept.

To avoid redundancy in the storage of points, the input points of a segment Delaunay graph
are stored in a container, and the various types of sites (input points and segments,
points of intersection, subsegments with one or two points of intersection as endpoints)
only store handles to the points in the container.

See section \ref Segment_Delaunay_graph_2StronglyIntersecting for more information.

\tparam Gt must be a model of the concept `SegmentDelaunayGraphTraits_2`.

\cgalModels{SegmentDelaunayGraphStorageTraits_2}

\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`
*/
template<class Gt>
class Segment_Delaunay_graph_storage_traits_2
{
public:
  /*!
  */
  typedef unspecified_type Storage_site_2;
};

} //namespace CGAL
