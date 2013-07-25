/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeletonFace_2` describes the requirements for the face type of the 
`StraightSkeleton_2` concept. It is a refinement of the `HalfedgeDSFace` concept 
with support for storage of the incident halfedge. 

\cgalRefines `HalfedgeDSFace` 

\cgalHasModel CGAL::Straight_skeleton_face_base_2

\sa `StraightSkeleton_2` 
\sa `CGAL::Straight_skeleton_face_base_2<Refs>` 
\sa `CGAL::Straight_skeleton_vertex_base_2<Refs,Point,FT>` 
\sa `CGAL::Straight_skeleton_halfedge_base_2<Refs>` 

*/

class StraightSkeletonFace_2 {
public:

/// \name Creation 
/// @{

/*!
Default constructor 
*/ 
StraightSkeletonFace_2(); 

/*!
Constructs a face with ID number `id` 
*/ 
StraightSkeletonFace_2(int id); 

/// @} 

/// \name Access Functions 
/// @{

/*!
The ID of the face. 
*/ 
int id() const; 

/// @}

}; /* end StraightSkeletonFace_2 */
