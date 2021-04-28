/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines TriangulationOnSphereTraits_2

The concept `DelaunayTriangulationOnSphereTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Delaunay_triangulation_on_sphere_2<Traits, Tds>`.

To the requirements listed within the concept `TriangulationOnSphereTraits_2`,
this concept adds types and functors requirements related to build the dual on the sphere.

\cgalHasModel `CGAL::Delaunay_triangulation_on_sphere_traits_2`
\cgalHasModel `CGAL::Projection_on_sphere_traits_3`
*/
class DelaunayTriangulationOnSphereTraits_2
{
public:
  /// \name Constructions
  ///
  /// @{

  /// Construction object type. Must provide the operator:
  ///
  /// `Point_on_sphere_2 operator()(Point_on_sphere_2 p, Point_on_sphere_2 q, Point_on_sphere_2 r)`
  ///
  /// which returns the intersection of the dual of the face defined by the three points `p`, `q`, and `r`,
  /// and the sphere, on the positive side of the plane defined by `p`, `q`, and `r`. The dual of the face
  /// is the line orthogonal to the face, passing through the center of the smallest circumscribing sphere
  /// of the face.
  ///
  /// \pre the center of the sphere is on the negative side of the plane defined by `p`, `q`, and `r`.
  ///
  /// \note This type is only required for the computation of dual objects (Voronoi vertices and edges)
  /// and a dummy type can be used otherwise.
  typedef unspecified_type Construct_circumcenter_on_sphere_2;

  /// Construction object type. Must provide the operator:
  ///
  /// `Point_3 operator()(Point_3 p, Point_3 q, Point_3 r)`
  ///
  /// which returns the center of the circle circumscribed to the face with vertices `p`, `q`, and `r`.
  ///
  /// \note This type is only required for the computation of dual objects (Voronoi vertices and edges)
  /// and a dummy type can be used otherwise.
  typedef unspecified_type Construct_circumcenter_3;

  /// @}

public:
  /// \name Operations
  ///
  /// The following functions give access to the predicate and constructor objects.
  ///
  /// @{

  ///
  Construct_circumcenter_on_sphere_2 construct_circumcenter_on_sphere_2_object();

  ///
  Construct_circumcenter_3 construct_circumcenter_3_object();

  /// @}
};
