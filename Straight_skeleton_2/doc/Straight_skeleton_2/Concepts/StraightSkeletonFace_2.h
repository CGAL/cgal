/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

\cgalRefines{HalfedgeDSFace}

The concept `StraightSkeletonFace_2` describes the requirements for the face type of the
`StraightSkeleton_2` concept. It is a refinement of the `HalfedgeDSFace` concept
with support for storage of the incident halfedge.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Straight_skeleton_face_base_2}
\cgalHasModelsEnd

\sa `StraightSkeletonVertex_2`
\sa `StraightSkeletonHalfedge_2`
*/
class StraightSkeletonFace_2 {
public:

/// \name Creation
/// @{

/*!
%Default constructor
*/
StraightSkeletonFace_2();

/*!
constructs a face with ID number `id`
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
