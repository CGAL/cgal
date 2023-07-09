
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_hierarchy_vertex_base_2` provides a model for the
`SegmentDelaunayGraphHierarchyVertexBase_2` concept, which is the
vertex base required by the
`Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>` class.

\tparam Vbb must be a model of the `SegmentDelaunayGraphVertexBase_2` concept.

\cgalModels{SegmentDelaunayGraphHierarchyVertexBase_2}

\sa `SegmentDelaunayGraphVertexBase_2`
\sa `SegmentDelaunayGraphHierarchyVertexBase_2`
\sa `SegmentDelaunayGraphDataStructure_2`
\sa `CGAL::Segment_Delaunay_graph_vertex_base_2<St,Vb>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>`

*/
template< typename Vbb >
class Segment_Delaunay_graph_hierarchy_vertex_base_2 : public Vbb {
public:

}; /* end Segment_Delaunay_graph_hierarchy_vertex_base_2 */
} /* end namespace CGAL */
