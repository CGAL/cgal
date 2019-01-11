
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraphLinf2Ref

The class `Segment_Delaunay_graph_Linf_2`
represents the segment Delaunay graph under the \f$ L_{\infty} \f$ metric (which is
the dual graph of the 2D segment Voronoi diagram under the \f$ L_{\infty} \f$
metric).
This class is derived from class `CGAL::Segment_Delaunay_graph_2<Gt,DS>`.
Currently, deletion of sites is not supported.

\tparam Gt is the geometric traits
which must be a model of `SegmentDelaunayGraphLinfTraits_2`

\tparam DS is the segment Delaunay graph data structure
which must be a model of `SegmentDelaunayGraphDataStructure_2`.
`DS` defaults to
`CGAL::Triangulation_data_structure_2< CGAL::Segment_Delaunay_graph_vertex_base_2<Gt>, CGAL::Triangulation_face_base_2<Gt> >`.

\cgalModels `DelaunayGraph_2`

\sa `SegmentDelaunayGraphLinfTraits_2`
\sa `SegmentDelaunayGraphDataStructure_2`
\sa `CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename Gt, typename DS >
class Segment_Delaunay_graph_Linf_2 {
}; /* end Segment_Delaunay_graph_Linf_2 */
} /* end namespace CGAL */
