/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines TriangulationOnSphereTraits_2
\cgalRefines SpatialSortingTraits_3

The concept `DelaunayTriangulationOnSphereTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Delaunay_triangulation_on_sphere_2<Traits, Tds>`.
This concept provides the types of the geometric primitives used in the
triangulation and some function object types for the required predicates on those primitives.

\cgalHasModel `CGAL::Delaunay_triangulation_sphere_traits_2<K>`
\cgalHasModel `CGAL::Geographical_coordinates_traits_2<K>`
\cgalHasModel `CGAL::Projection_sphere_traits_3<K>`

*/
class DelaunayTriangulationOnSphereTraits_2
{
  /*!
  The 3D segment type.
  */
  typedef unspecified_type Segment_3;

public:
  /// Construction object type. Must provide the operator
  ///
  /// `Point_on_sphere_2 operator()(Point_on_sphere_2 p, Point_on_sphere_2 q, Point_on_sphere_2 r)`
  ///
  /// which returns the intersection of the dual of the facet defined by the three points `p`, `q`, and `r`,
  /// and the sphere.
  ///
  /// \note This type is only required for the computation of dual objects
  /// and a dummy type can be used otherwise.
  typedef unspecified_type Construct_circumcenter_on_sphere_2;

  /// Construction object which must provide the following operator
  ///
  /// `Segment_3 operator()(Point_on_sphere_2 p, Point_on_sphere_2 q)`
  ///
  /// which introduces an arc of great circle, with source `p` and target `q`.
  ///
  /// \note This type is only required for the computation of dual objects
  /// and a dummy type can be used otherwise.
  typedef unspecified_type Construct_segment_3;

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
  Construct_segment_3 construct_segment_3_object();

  /// @}
};
