
namespace CGAL {

/*!
\ingroup PkgSDGLinf

The class `Segment_Delaunay_graph_Linf_2`
represents the segment Delaunay graph under the L_{infinity} metric (which is
the dual graph of the 2D segment Voronoi diagram under the L_{infinity}
metric).
Currently it supports only insertions of sites.
It is templated by two template arguments `Gt`, which
must be a model of `SegmentDelaunayGraphLinfTraits_2`
and `DS`,
which must be a model of `SegmentDelaunayGraphDataStructure_2`.
The second template argument defaults to
`CGAL::Triangulation_data_structure_2< CGAL::Segment_Delaunay_graph_vertex_base_2<Gt>, CGAL::Triangulation_face_base_2<Gt> >`.
This class is derived from class
`CGAL::Segment_Delaunay_graph_2<Gt,DS>`.

\cgalModels `DelaunayGraph_2`

\sa `DelaunayGraph_2`
\sa `SegmentDelaunayGraphLinfTraits_2`
\sa `SegmentDelaunayGraphDataStructure_2`
\sa `SegmentDelaunayGraphVertexBase_2`
\sa `TriangulationFaceBase_2`
\sa `CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Segment_Delaunay_graph_vertex_base_2<Gt,SSTag>`
\sa `CGAL::Triangulation_face_base_2<Gt>`

*/
template< typename Gt, typename DS >
class Segment_Delaunay_graph_Linf_2 {
}; /* end Segment_Delaunay_graph_Linf_2 */
} /* end namespace CGAL */
