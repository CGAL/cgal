/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeletonItemsConverter_2` describes the requirements for items converter passed
as the third template argument to the class `Straight_skeleton_converter_2<SrcSs,TgtSs,ItemsConverter>`.
It converts the `HDS` items from one type of straight skeleton to another.

\cgalHasModel `CGAL::Straight_skeleton_items_converter_2`

\sa `CGAL::Straight_skeleton_converter_2<SrcSs,TgtSs,ItemsCvt>`
*/
class StraightSkeletonItemsConverter_2 {
public:

/// \name Types
/// @{

/*!
A constant handle to a model of the `StraightSkeletonVertex_2` concept
used to represent the vertices of the straight skeleton of the source type.
*/
typedef unspecified_type Source_vertex_const_handle;

/*!
A constant handle to model of the `StraightSkeletonHalfedge_2` concept
used to represent the halfedges of the straight skeleton of the source type.
*/
typedef unspecified_type Source_halfedge_const_handle;

/*!
A constant handle to model of the `StraightSkeletonFace_2` concept
used to represent the faces of the straight skeleton of the source type.
*/
typedef unspecified_type Source_face_const_handle;

/*!
A model of the `StraightSkeletonVertex_2` concept used to represent the vertices
of the straight skeleton of the target type.
*/
typedef unspecified_type Target_vertex;

/*!
A model of the `StraightSkeletonHalfedge_2` concept used to represent the halfedges
of the straight skeleton of the target type.
*/
typedef unspecified_type Target_halfedge;

/*!
Any model of the `StraightSkeletonFace_2` concept used to represent the faces
of the straight skeleton of the target type.
*/
typedef unspecified_type Target_face;

/// @}

/// \name Operations
/// @{

/*!
returns a new vertex with the same data as `v` converted to the corresponding target types.
*/
Target_vertex operator()( Source_vertex_const_handle v) const;

/*!
returns a new halfedge with the same data as `h` converted to the corresponding target types.
*/
Target_halfedge operator()( Source_halfedge_const_handle h) const;

/*!
returns a new face with the same data as `f` converted to the corresponding target types.
*/
Target_face operator()( Source_face_const_handle f) const;

/// @}

}; /* end StraightSkeletonItemsConverter_2 */
