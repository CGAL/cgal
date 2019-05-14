
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2

The class `Segment_Delaunay_graph_vertex_base_2` provides a model for the 
`SegmentDelaunayGraphVertexBase_2` concept which is the vertex 
base required by the `SegmentDelaunayGraphDataStructure_2` 
concept. 

\tparam Gt must be a model of the concept `SegmentDelaunayGraphTraits_2`. 
\tparam SSTag indicates whether 
or not to use the simple storage site that does not support 
intersecting segments, or the full storage site, that supports 
intersecting segments. The possible values are `Tag_true` 
and `Tag_false`. `Tag_true` indicates that the 
full storage site is to be used, whereas `Tag_false` 
indicates that the simple storage site is to be used. 

\cgalModels `SegmentDelaunayGraphVertexBase_2`

\sa `SegmentDelaunayGraphVertexBase_2` 
\sa `SegmentDelaunayGraphDataStructure_2` 
\sa `SegmentDelaunayGraphTraits_2` 
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>` 
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>` 

*/
template< typename Gt, typename SSTag >
class Segment_Delaunay_graph_vertex_base_2 {
public:

/// @}

}; /* end Segment_Delaunay_graph_vertex_base_2 */
} /* end namespace CGAL */
