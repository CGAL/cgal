/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines `TriangulationDSFaceBase_2`

The concept `TriangulationOnSphereFaceBase_2` describes the requirements for a face class
of a triangulation data structure that is itself plugged into a triangulation on the sphere.

The data structure concept `TriangulationDataStructure_2` was primarily designed
to serve as a data structure for the 2D triangulation classes of \cgal, which are triangulations
embedded in the 2D Euclidean plane.
However, its genericity makes it usable for any orientable triangulated surface without boundary
regardless of the dimensionality of the space the triangulation is embedded in, and thus
it is a valid data structure for the triangulations on the sphere of this package.

A departing feature is that if the vertices of the triangulation are all lying on one hemisphere,
then the triangulation is not an orientable triangulated surface without boundary.
In this case, fictitious faces are added to the triangulation, called <i>ghost faces</i>,
such that the triangulation is a topological sphere.

\cgalHasModel `CGAL::Triangulation_on_sphere_face_base_2<Traits, Fb>`

\sa `TriangulationDataStructure_2`
\sa `TriangulationOnSphereVertexBase_2`
*/
class TriangulationOnSphereFaceBase_2
{
public:
  /// provides read access to a Boolean used to indicate if the face is a ghost face.
  bool is_ghost();

  /// provides write access to a Boolean used to indicate if the face is a ghost face.
  void set_ghost(const bool b);
};
