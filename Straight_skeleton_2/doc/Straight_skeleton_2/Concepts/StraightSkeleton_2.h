/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeleton_2` describes the requirements for the data structure 
used to represent a straight skeleton. It refines the concept 
`HalfedgeDS` and adds additional requirements on the nested types 
`Vertex`, `Halfedge`, and `Face` of the halfedge data structure. 

\cgalRefines `HalfedgeDS` 

\cgalHasModel CGAL::Straight_skeleton_2

\attention This concept explicitly protects all the modifying
operations of the base `HalfedgeDS` concept. Only the algorithm
classes, or clients explicitly bypassing the protection mechanism, can
modify a straight skeleton.

*/

class StraightSkeleton_2 {
public:

/// \name Types 
/// @{

/*!
A model of the `StraightSkeletonVertex_2` concept used to represent the vertices of the straight skeleton 
*/ 
typedef unspecified_type Vertex; 

/*!
A model of the `StraightSkeletonHalfedge_2` concept used to represent the halfedges of the straight skeleton 
*/ 
typedef unspecified_type Halfedge; 

/*!
Any model of the `StraightSkeletonFace_2` concept used to represent the faces of the straight skeleton 
*/ 
typedef unspecified_type Face; 

/// @}

}; /* end StraightSkeleton_2 */
