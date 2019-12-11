/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines SpatialSortingTraits_2

The concept `DelaunayTriangulationOnSphereTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Delaunay_triangulation_on_sphere_2<Traits, Tds>`.
This concept provides the types of the geometric primitives used in the
triangulation and some function object types for the required predicates on those primitives.

In particular, the traits class is expected to contain information about the sphere (center and radius)
on which the points of the triangulation live, as well as provide means to change this information.

@todo explain too close points

\cgalHasModel `CGAL::Delaunay_triangulation_sphere_traits_2<K>`
\cgalHasModel `CGAL::Projection_sphere_traits_3<K>`

*/
class DelaunayTriangulationOnSphereTraits_2
{
public:
  /// \name Types
  /// @{

  /*!
  The number type; must be a model of `FieldNumberType`.
  */
  typedef unspecified_type FT;

  /*!
  The point type of the triangulation, representing a point on the sphere.
  */
  typedef unspecified_type Point_on_sphere;

  /*!
  The 3D point type.
  */
  typedef unspecified_type Point_3;

  /*!
  The 3D segment type.
  */
  typedef unspecified_type Segment_3;

public:
  // Predicates and constructions

  /// Predicate object which must provide the following operator
  typedef unspecified_type Compare_xyz_3;

  /// Predicate object
  typedef unspecified_type Equal_on_sphere_2;

  /// Predicate object
  typedef unspecified_type Inside_cone_2;

  /// Predicate object
  typedef unspecified_type Orientation_on_sphere_2;

  /// Predicate object
  typedef unspecified_type Orientation_3;

  /// Predicate object
  typedef unspecified_type Power_test_2;

public:
  // Only needed for dual()

  /// Construction object which must provide the following operator
  /// only required for dual() functions
  typedef unspecified_type Construct_circumcenter_on_sphere_2;

  /// Construction object which must provide the following operator
  /// only required for dual() functions
  typedef unspecified_type Construct_segment_3;

public:
  // Really specific to ToS

  /// Returns whether the point is on the sphere or not
  /// A point not on the sphere is not inserted
  bool is_on_sphere(Point_on_sphere); // @todo a bit awkward

  /// Returns the center of the sphere.
  const Point_3& center();

  /// Sets the center of the sphere.
  void set_center(Point_3);

  /// Returns the radius of the sphere.
  FT radius();

  /// Sets the radius of the sphere.
  void set_radius(FT radius);

public:
  /// Returns whether two points are too close
  bool are_points_too_close();
};
