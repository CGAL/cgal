namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Auxiliary

The class `Straight_skeleton_halfedge_base_2` provides a model for the 
`StraightSkeletonHalfedge_2` concept which is the halfedge 
type required by the `StraightSkeleton_2` 
concept.
\tparam Refs a model of the `StraightSkeleton_2` concept (the halfedge container)

This class can be used as a base class allowing users of the straight skeleton data structure to decorate a halfedge with additional data. The concrete halfedge class must be given in the `HalfedgeDSItems` template parameter of the instantiation of the `HalfedgeDS_default` class used as the model for the `Straight_skeleton_2` concept. 

\cgalModels `StraightSkeletonHalfedge_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `StraightSkeletonHalfedge_2` 
\sa `StraightSkeletonVertex_2` 
\sa `StraightSkeleton_2` 
\sa `Straight_skeleton_vertex_base_2`

*/
template< typename Refs >
class Straight_skeleton_halfedge_base_2 {
public:

/// @}

}; /* end Straight_skeleton_halfedge_base_2 */
} /* end namespace CGAL */
