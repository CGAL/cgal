
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraphLinf2Ref

The class `Segment_Delaunay_graph_Linf_traits_2` provides a model for the
`SegmentDelaunayGraphLinfTraits_2` concept.
The class is derived from class
`CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`.
We refer to the documentation of the base class
`CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
for an explanation of the template parameters.
Our experience has shown that the
`CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
class is more efficient for all
reasonable and valid values of the template parameters
and should be preferred over
`Segment_Delaunay_graph_Linf_traits_2`.

\cgalModels `SegmentDelaunayGraphLinfTraits_2`

\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename K, typename MTag >
struct Segment_Delaunay_graph_Linf_traits_2 {
}; /* end Segment_Delaunay_graph_Linf_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraphLinf2Ref

The class `Segment_Delaunay_graph_Linf_traits_without_intersections_2`
provides a model for the
`SegmentDelaunayGraphLinfTraits_2` concept.
The class is derived from class
`CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`.
We refer to the documentation of the base class
`CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
for an explanation of the template parameters.
Our experience has shown that the
`CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`
class is more efficient for all
reasonable and valid values of the template parameters
and should be preferred over
`Segment_Delaunay_graph_Linf_traits_without_intersections_2`.

\cgalModels `SegmentDelaunayGraphLinfTraits_2`

\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename K, typename MTag >
struct Segment_Delaunay_graph_Linf_traits_without_intersections_2 {
}; /* end Segment_Delaunay_graph_Linf_traits_without_intersections_2 */
} /* end namespace CGAL */
