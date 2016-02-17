namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Auxiliary

The class `Straight_skeleton_face_base_2` provides a model for the 
`StraightSkeletonFace_2` concept which is the face 
type required by the `StraightSkeleton_2` 
concept.
\tparam Refs a model of the `StraightSkeleton_2` concept (the face container)

This class can be used as a base class allowing users of the straight skeleton data structure to decorate a face with additional data. The concrete face class must be given in the `HalfedgeDSItems` template parameter of the instantiation of the `HalfedgeDS_default` class used as the model for the `Straight_skeleton_2` concept. 

\cgalModels `StraightSkeletonFace_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `StraightSkeletonFace_2` 
\sa `StraightSkeletonVertex_2` 
\sa `StraightSkeletonHalfedge_2` 
\sa `StraightSkeleton_2` 
\sa `Straight_skeleton_vertex_base_2` 
\sa `Straight_skeleton_halfedge_base_2` 

*/
template< typename Refs >
class Straight_skeleton_face_base_2 {
public:

/// @}

}; /* end Straight_skeleton_face_base_2 */
} /* end namespace CGAL */
