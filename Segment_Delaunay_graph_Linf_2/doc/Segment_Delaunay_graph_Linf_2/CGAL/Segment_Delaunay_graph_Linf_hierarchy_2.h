
namespace CGAL {

/*!
\ingroup PkgSDGLinf

This class is equivalent to the `Segment_Delaunay_graph_hierarchy_2` class,
but it uses `Segment_Delaunay_graph_Linf_2<Gt,DS>`
at every level of the hierarchy
(instead of `Segment_Delaunay_graph_2<Gt,DS>` at each level of
`Segment_Delaunay_graph_hierarchy_2`).
For more details related to the template parameters, see the documentation
of `Segment_Delaunay_graph_hierarchy_2`.

The `Segment_Delaunay_graph_Linf_hierarchy_2` class derives publicly from
the `Segment_Delaunay_graph_Linf_2<Gt,DS>` class.
The interface is the same as its base class.

\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `SegmentDelaunayGraphDataStructure_2`
\sa `SegmentDelaunayGraphLinfTraits_2`
\sa `SegmentDelaunayGraphHierarchyVertexBase_2`
\sa `CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<Vbb>`

*/
template< typename Gt, typename STag, typename DS >
class Segment_Delaunay_graph_Linf_hierarchy_2 :
  public CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS> {
}; /* end Segment_Delaunay_graph_Linf_hierarchy_2 */
} /* end namespace CGAL */
