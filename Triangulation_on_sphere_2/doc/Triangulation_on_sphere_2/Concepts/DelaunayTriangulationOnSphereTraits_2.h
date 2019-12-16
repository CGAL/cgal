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
  /// `Point_3 operator()(Point_on_sphere p, Point_on_sphere q, Point_on_sphere r)`
  ///
  /// which returns the center of the smallest ball passing through the three points `p`, `q`, and `r`.
  ///
  /// @todo should this actually be something like `Construct_circumcenter_on_sphere_2`
  /// and project the circumcenter on the sphere rather than just computing the center of the smallest
  /// circumscribing ball of the triangle?
  ///
  /// \note This type is only required for the computation of dual objects
  /// and a dummy type can be used otherwise.
  typedef unspecified_type Construct_circumcenter_3;

  /// Construction object which must provide the following operator
  ///
  /// `Point_3 operator()(Point_on_sphere p)`
  ///
  /// which constructs a 3D point from the point on the sphere.
  ///
  /// \note This operator should return a const reference if it makes sense, to avoid superfluous copies.
  typedef unspecified_type Construct_point_3;

  /// Construction object which must provide the following operator
  ///
  /// `Segment_3 operator()(Point_3 p, Point_3 q)`
  ///
  /// which introduces a segment with source `p` and target `q`.
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
  Construct_circumcenter_3 construct_circumcenter_3_object();

  ///
  Construct_segment_3 construct_segment_3_object();

  /// @}
};
