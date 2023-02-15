
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

\cgalModels `SegmentDelaunayGraphVertexBase_2`

The class `Segment_Delaunay_graph_vertex_base_2` provides a model for the
`SegmentDelaunayGraphVertexBase_2` concept which is the vertex
base required by the `SegmentDelaunayGraphDataStructure_2`
concept.

\tparam St must be a model of the concept `SegmentDelaunayGraphStorageTraits_2`.
           This type must be template parameter used for `CGAL::Segment_Delaunay_graph_2`.

\tparam Vb is a vertex base class from which `Segment_Delaunay_graph_vertex_base_2` derives.
           It must be a model of the `TriangulationDSVertexBase_2` concept.
           It has the default value `CGAL::Triangulation_ds_vertex_base_2<>`.

\sa `CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<Vb>`
*/
template< typename St, typename Vb >
class Segment_Delaunay_graph_vertex_base_2 {
public:

}; /* end Segment_Delaunay_graph_vertex_base_2 */
} /* end namespace CGAL */
