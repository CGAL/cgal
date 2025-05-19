namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Auxiliary

\cgalModels{StraightSkeletonVertex_2}

The class `Straight_skeleton_vertex_base_2` provides a model for the `StraightSkeletonVertex_2` concept
which is the vertex type required by the `StraightSkeleton_2` concept.

\tparam Refs must be a model of `StraightSkeleton_2`
\tparam Point a Point type
\tparam FT must be a model of the `FieldWithSqrt`, which is the numeric type used to represent
           the time of a vertex (a Euclidean distance).

This class can be used as a base class allowing users of the straight skeleton data structure to decorate
a vertex with additional data. The concrete vertex class must be given in the `HalfedgeDSItems`
template parameter of the instantiation of the `HalfedgeDS_default` class used as the model
for the `CGAL::Straight_skeleton_2` concept.

\sa `StraightSkeletonFace_2`
\sa `StraightSkeletonVertex_2`
\sa `StraightSkeletonHalfedge_2`
\sa `StraightSkeleton_2`
\sa `CGAL::Straight_skeleton_halfedge_base_2`
\sa `CGAL::Straight_skeleton_face_base_2`
*/
template< typename Refs, typename Point, typename FT >
class Straight_skeleton_vertex_base_2 {
public:

}; /* end Straight_skeleton_vertex_base_2 */

} /* end namespace CGAL */
