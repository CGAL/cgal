
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_traits_2` provides a model for the
`SegmentDelaunayGraphTraits_2` concept.

\tparam K  must be a model of the `Kernel` concept.

\tparam MTag corresponds to how predicates are evaluated. There are two
possible values for `MTag`, namely `Field_with_sqrt_tag` and
`Field_tag`. The first one must be used when the number type
used in the representation supports the exact evaluation of signs of
expressions involving all four basic operations and square roots,
whereas the second one requires the exact evaluation of signs of
field-type expressions, i.e., expressions involving additions,
subtractions, multiplications and divisions. The default value for
`MTag` is `Field_tag`.
The way the predicates are evaluated is discussed in
\cgalCite{b-ecvdl-96} and \cgalCite{cgal:k-reisv-04} (the geometric filtering
part).

\cgalModels `SegmentDelaunayGraphTraits_2`

\sa `Kernel`
\sa `SegmentDelaunayGraphTraits_2`
\sa `CGAL::Field_tag`
\sa `CGAL::Field_with_sqrt_tag`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,St,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename K, typename MTag >
struct Segment_Delaunay_graph_traits_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_true Intersections_tag;

/// @}

/// \name Additional Types
/// The `Segment_Delaunay_graph_traits_2` class introduces a few
/// additional types with respect to the
/// `SegmentDelaunayGraphTraits_2` concept. These are:
/// @{

/*!
A typedef for the template parameter
`K`.
*/
typedef K Kernel;

/*!
A typedef for the template
parameter `MTag`.
*/
typedef MTag Method_tag;

/// @}

}; /* end Segment_Delaunay_graph_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_traits_without_intersections_2` provides a model for the
`SegmentDelaunayGraphTraits_2` concept.

\tparam K  must be a model of the `Kernel`
concept.

\tparam MTag corresponds to how predicates
are evaluated. There are two possible values for `MTag`, namely
`Field_with_sqrt_tag` and `Euclidean_ring_tag`. The first one
must be used when the number type used in the representation supports
the exact evaluation of signs of expressions involving all four basic
operations and square roots, whereas the second one requires the exact
evaluation of signs of ring-type expressions, i.e., expressions
involving only additions, subtractions and multiplications. The
default value for `MTag` is `Euclidean_ring_tag`.
The way the predicates are evaluated is discussed in
\cgalCite{b-ecvdl-96} and \cgalCite{cgal:k-reisv-04} (the geometric filtering
part).

\cgalModels `SegmentDelaunayGraphTraits_2`

\sa `Kernel`
\sa `SegmentDelaunayGraphTraits_2`
\sa `CGAL::Euclidean_ring_tag`
\sa `CGAL::Field_with_sqrt_tag`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,St,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,STag,DS>`
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>`

*/
template< typename K, typename MTag >
struct Segment_Delaunay_graph_traits_without_intersections_2 {

/// \name Types
/// @{

/*!

*/
typedef CGAL::Tag_false Intersections_tag;

/// @}

/// \name Additional Types
/// The `Segment_Delaunay_graph_traits_without_intersections_2` class
/// introduces a few additional types with respect to the
/// `SegmentDelaunayGraphTraits_2` concept. These are:
/// @{

/*!
A typedef for the template parameter
`K`.
*/
typedef K Kernel;

/*!
A typedef for the template
parameter `MTag`.
*/
typedef MTag Method_tag;

/// @}

}; /* end Segment_Delaunay_graph_traits_without_intersections_2 */
} /* end namespace CGAL */
