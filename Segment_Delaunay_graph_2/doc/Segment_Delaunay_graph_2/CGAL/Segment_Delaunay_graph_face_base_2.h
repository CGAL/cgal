
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

\cgalModels{SegmentDelaunayGraphFaceBase_2}

The class `Segment_Delaunay_graph_face_base_2` provides a model for the
`SegmentDelaunayGraphFaceBase_2` concept which is the face
base required by the `SegmentDelaunayGraphDataStructure_2`
concept.

\tparam Gt must be a model of the concept `SegmentDelaunayGraphTraits_2`.
           This type must be identical to the template parameter used for `CGAL::Segment_Delaunay_graph_2`.

\tparam Fb is a face base class from which `Segment_Delaunay_graph_face_base_2` derives.
           It must be a model of the `TriangulationFaceBase_2` concept.
           It has the default value `CGAL::Triangulation_face_base_2<Gt>`.
*/
template< typename Gt, typename Fb >
class Segment_Delaunay_graph_face_base_2 {
public:

}; /* end Segment_Delaunay_graph_face_base_2 */

} /* end namespace CGAL */
