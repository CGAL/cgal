/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines `TriangulationDSFaceBase_2`

The concept `TriangulationOnSphereFaceBase_2` describes the requirements for a face class
of a triangulation data structure that is itself plugged into a Delaunay triangulation
of the sphere.

The data structure concept `TriangulationDataStructure_2` was primarily designed
to serve as a data structure for the 2D triangulation classes of \cgal, which are triangulations
embedded in the 2D Euclidean plane.
However its genericy makes it usable for any orientable triangulated surface without boundary,
whatever may be the dimensionality of the space the triangulation is embedded in, and thus
it is a valid data structure for the triangulations on the sphere of this package.

A detail is that if the vertices of the triangulation are all lying on one hemisphere,
then the triangulation is not an orientable triangulated surface without boundary.
In this case, the triangulations on sphere add fictitious faces, called <i>ghost faces</i>,
such that the triangulation is topologicial sphere.

Another Boolean is required to mark whether a face is in conflict or not, for internal machinery.

\cgalHasModel `CGAL::Triangulation_on_sphere_face_base_2<Traits, Fb>`

\sa `TriangulationDataStructure_2`
\sa `TriangulationOnSphereVertexBase_2`
*/
class TriangulationOnSphereFaceBase_2
{
public:
  /// Provides access to the ghost Boolean so that it may be read or modified.
  bool& ghost();

  /// Provides access to the conflict Boolean so that it may read or modified.
  bool& in_conflict();
};
