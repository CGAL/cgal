namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Auxiliary

The class `Straight_skeleton_vertex_base_2` provides a model for the 
`StraightSkeletonVertex_2` concept which is the vertex 
type required by the `StraightSkeleton_2` 
concept. 
\tparam Refs a model of the `StraightSkeleton_2` concept (the vertex container)
\tparam Point a Point type
\tparam FT a model of the `FieldWithSqrt` concept, which is the numeric type used to represent the time of a vertex (a Euclidean distance). 

This class can be used as a base class allowing users of the straight skeleton data structure to decorate a vertex with additional data. The concrete vertex class must be given in the `HalfedgeDSItems` template parameter of the instantiation of the `HalfedgeDS_default` class used as the model for the `Straight_skeleton_2` concept. 

\cgalModels `StraightSkeletonVertex_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `StraightSkeletonVertex_2` 
\sa `StraightSkeletonHalfedge_2` 
\sa `StraightSkeleton_2` 
\sa `CGAL::Straight_skeleton_halfedge_base_2`

*/
template< typename Refs, typename Point, typename FT >
class Straight_skeleton_vertex_base_2 {
public:

/// @}

}; /* end Straight_skeleton_vertex_base_2 */
} /* end namespace CGAL */
