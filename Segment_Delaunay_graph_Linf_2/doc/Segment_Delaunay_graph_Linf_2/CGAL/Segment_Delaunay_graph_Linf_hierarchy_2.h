
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraphLinf2Ref

This class is equivalent to the `Segment_Delaunay_graph_hierarchy_2` class,
but it uses `Segment_Delaunay_graph_Linf_2<Gt,DS>`
at every level of the hierarchy
(instead of `Segment_Delaunay_graph_2<Gt,DS>` at each level of
`Segment_Delaunay_graph_hierarchy_2`).

The `Segment_Delaunay_graph_Linf_hierarchy_2` class derives publicly from
the `Segment_Delaunay_graph_Linf_2<Gt,DS>` class and has the same
interface.
For more details related to the template parameters, see the documentation
of `Segment_Delaunay_graph_hierarchy_2`.

The class `Segment_Delaunay_graph_Linf_hierarchy_2`
should be preferred when the size of the input is large
and the input sites are inserted one by one (e.g., when the
whole input is not known in advance):
Experiments suggest that the hierarchy is faster than the plain
segment Delaunay graph when the size of the input exceeds 1000 sites.
If the whole input is known in advance, another option is to use
the plain segment Delaunay graph but only after
spatial sorting of the sites.
See also the examples in the User Manual for different uses of the hierarchy
and spatial sorting.

\sa `SegmentDelaunayGraphLinfTraits_2`
\sa `CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename Gt, typename STag, typename DS >
class Segment_Delaunay_graph_Linf_hierarchy_2 :
  public CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS> {
}; /* end Segment_Delaunay_graph_Linf_hierarchy_2 */
} /* end namespace CGAL */
