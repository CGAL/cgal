/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

The concept `TriangulationOnSphereFaceBase_2` describes the requirements for the base face class
of a triangulation data structure that is itself plugged into a Delaunay triangulation
of the sphere.

This concept refines the concept `TriangulationDSFaceBase_2` by adding two functions
related to so-called <i>ghost faces</i>. Ghost faces are faces of the triangulation data structure
that ...
@todo FINISH WRITING THIS
@todo DO I NEED SET_CONFLICT_FLAG ?

\cgalRefines `TriangulationDSFaceBase_2`

\cgalHasModel `CGAL::Triangulation_on_sphere_face_base_2<Traits, Fb>`

\sa `TriangulationDataStructure_2`
\sa `TriangulationOnSphereVertexBase_2`
*/
class TriangulationOnSphereFaceBase_2
{
public:

  /// Returns whether the face is a <i>ghost</i> face, or not.
  bool is_ghost() const;

  /// Provides access to the ghost Boolean so that it may be modified.
  bool& ghost();
};
