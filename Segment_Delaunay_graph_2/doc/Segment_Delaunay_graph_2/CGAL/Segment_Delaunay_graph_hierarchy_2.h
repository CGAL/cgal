
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

We provide an alternative to the class
`Segment_Delaunay_graph_2<Gt,St,DS>` for the incremental
construction of the segment Delaunay graph.
The `Segment_Delaunay_graph_hierarchy_2` class
maintains a hierarchy of Delaunay graphs. There are two possibilities
as to how this hierarchy is constructed.

In the first case the bottom-most level of the hierarchy contains the
full segment Delaunay graph. The upper levels of the hierarchy
contain only points that are either point sites or endpoints of
segment sites in the bottom-most Delaunay graph.
A point that is in level \f$ i\f$ (either as an individdual point or as the
endpoint of a segment), is inserted in level \f$ i+1\f$ with probability
\f$ 1/\alpha\f$ where \f$ \alpha>1\f$ is some constant.
In the second case the upper levels of the hierarchy contains not only
points but also segments. A site that is in level \f$ i\f$, is in level
\f$ i+1\f$ with probability \f$ 1/\beta\f$ where \f$ \beta > 1\f$ is some constant.

The difference between the `Segment_Delaunay_graph_2<Gt,St,DS>`
class and the `Segment_Delaunay_graph_hierarchy_2` class (both versions of it) is on how the
nearest neighbor location is done. Given a point \f$ p\f$ the location is
done as follows: at the top most level we find the nearest neighbor of
\f$ p\f$ as in the `Segment_Delaunay_graph_2<Gt,St,DS>` class. At
every subsequent level \f$ i\f$ we use the nearest neighbor found at level
\f$ i+1\f$ to find the nearest neighbor at level \f$ i\f$. This is a variant of
the corresponding hierarchy for points found in \cgalCite{cgal:d-dh-02}. The
details are described in \cgalCite{cgal:k-reisv-04}.

The class has four template parameters. The first and fourth
have essentially the same semantics as in the `Segment_Delaunay_graph_2<Gt,St,DS>` class.

\tparam Gt must be a model of the
`SegmentDelaunayGraphTraits_2` concept.

\tparam St must be a model of `SegmentDelaunayGraphStorageTraits_2`.
        By default, the storage traits is instantiated by `Segment_Delaunay_graph_storage_traits_2<Gt>`.

\tparam STag The third template
parameter controls whether or not segments are added in the upper
levels of the hierarchy. It's possible values are `Tag_true`
and `Tag_false`. If it is set to `Tag_true`,
segments are also inserted in the upper levels of the hierarchy. The
value `Tag_false` indicates that only points are to be
inserted in the upper levels of the hierarchy. The default value for
the third template parameter is `Tag_false`.

\tparam DS must be a model of the
`SegmentDelaunayGraphDataStructure_2` concept. However, the
vertex base class that is to be used in the segment Delaunay graph
data structure must be a model of the
`SegmentDelaunayGraphHierarchyVertexBase_2`
concept. The fourth template parameter defaults to:
\code
CGAL::Triangulation_data_structure_2<
          CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<
            CGAL::Segment_Delaunay_graph_vertex_base_2<St>,
          CGAL::Segment_Delaunay_graph_face_base_2<Gt> >
\endcode

The `Segment_Delaunay_graph_hierarchy_2` class derives publicly from the
`Segment_Delaunay_graph_2<Gt,St,DS>` class. The interface is
the same with its base class. In the sequel only additional types
and methods defined are documented.

\cgalModels{DefaultConstructible,CopyConstructible,Assignable}

\sa `CGAL::Segment_Delaunay_graph_2<Gt,St,DS>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<Vbb>`

*/
template< typename Gt, typename St, typename STag, typename DS >
class Segment_Delaunay_graph_hierarchy_2
  : public CGAL::Segment_Delaunay_graph_2<Gt,St,DS> {
public:

/// \name Types
/// `Segment_Delaunay_graph_hierarchy_2` introduces the following
/// types in addition to those introduced by its base class
/// `Segment_Delaunay_graph_2<Gt,St,DS>`.
/// @{

/*!
A type for the
`STag` template parameter.
*/
typedef STag Segments_in_hierarchy_tag;

/*!
A type for the base class.
*/
typedef CGAL::Segment_Delaunay_graph_2<Gt,St,DS> Base;

/// @}

/// \name Creation
/// In addition to the default and copy constructors, the following
/// constructors are defined:
/// @{

/*!
Creates a hierarchy of segment Delaunay graphs using
`gt` as geometric traits.
*/
Segment_Delaunay_graph_hierarchy_2(Gt gt=Gt());

/*!
Creates a segment Delaunay graph hierarchy using
`gt` as geometric traits and inserts all sites in the
range [`first`, `beyond`). `Input_iterator` must be a
model of `InputIterator`. The value type of `Input_iterator`
must be either `Point_2` or `Site_2`.
*/
template< class Input_iterator >
Segment_Delaunay_graph_hierarchy_2(Input_iterator first, Input_iterator beyond, Gt gt=Gt());

/// @}

}; /* end Segment_Delaunay_graph_hierarchy_2 */

/*!
Writes the current state of the segment Delaunay graph hierarchy to
an output stream. In particular, all sites in the diagram are
written to the stream (represented through appropriate input sites),
as well as the underlying combinatorial hierarchical data structure.
\relates Segment_Delaunay_graph_hierarchy_2
*/
std::ostream& operator<<(std::ostream& os, const Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>& svdh);

/*!
Reads the state of the segment Delaunay graph hierarchy from an
input stream.
\relates Segment_Delaunay_graph_hierarchy_2
*/
std::istream& operator>>(std::istream& is, const Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>& svdh);

} /* end namespace CGAL */
