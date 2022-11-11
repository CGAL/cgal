/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines SpatialSortingTraits_3

The concept `TriangulationOnSphereTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Triangulation_on_sphere_2<Traits, Tds>`.
This concept provides the types of the geometric primitives used in the
triangulation and the function object types for the required predicates on those primitives.

In particular, the traits class is expected to contain information about the sphere (center and radius)
on which the points of the triangulation lie, as well as to provide means to change this information.

\cgalHasModel `CGAL::Delaunay_triangulation_on_sphere_traits_2`
\cgalHasModel `CGAL::Projection_on_sphere_traits_3`

\sa `DelaunayTriangulationOnSphereTraits_2`
*/
class TriangulationOnSphereTraits_2
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
  typedef unspecified_type Point_on_sphere_2;

  /*!
  An arc of a great circle, used to represent a curved segment (a Voronoi or a Delaunay edge).
  */
  typedef unspecified_type Arc_on_sphere_2;

  /*!
  The 3D point type.
  */
  typedef unspecified_type Point_3;

  /*!
  The 3D segment type.
  */
  typedef unspecified_type Segment_3;

  /*!
  The 3D triangle type.
  */
  typedef unspecified_type Triangle_3;

  /// @}

public:
  /// \name Predicates
  /// @{

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere_2 p, Point_on_sphere_2 q, Point_on_sphere_2 r)`
  ///
  /// which returns `true` if `r` strictly lies on the shortest path between `p` and `q`
  /// on the great circle passing through `p` and `q`.
  ///
  /// If `p` and `q` are diametrically opposed, `true` is returned.
  typedef unspecified_type Collinear_are_strictly_ordered_on_great_circle_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere_2 p, Point_on_sphere_2 q)`
  ///
  /// which returns `true` or `false` based on a consistent, user-defined strict total order of points.
  typedef unspecified_type Compare_on_sphere_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere_2 p, Point_on_sphere_2 q)`
  ///
  /// which returns `true` if, and only if, `p` and `q` are equal.
  typedef unspecified_type Equal_on_sphere_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `Orientation operator()(Point_on_sphere_2 p, Point_on_sphere_2 q, Point_on_sphere_2 r)`
  ///
  /// which returns `CGAL::POSITIVE`, if `r` lies on the left hemisphere while walking
  /// the shortest path between `p` and `q` on the great circle defined by `p` and `q`;
  /// returns `CGAL::NEGATIVE` if `r` lies on the right hemisphere; and returns `CGAL::COPLANAR`
  /// if `r` lies on the great circle.
  ///
  /// If `p` and `q` are diametrically opposed, `CGAL::COPLANAR` is returned.
  typedef unspecified_type Orientation_on_sphere_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere_2 p, Point_on_sphere_2 q, Point_on_sphere_2 r, Point_on_sphere_2 s)`
  ///
  /// which returns `CGAL::POSITIVE`, if `s` lies on the positive side of the oriented plane `h`
  /// defined by `p`, `q`, and `r`; returns `CGAL::NEGATIVE` if `s` lies on the negative side of `h`;
  /// and returns `CGAL::COPLANAR` if `s` lies on `h`.
  typedef unspecified_type Side_of_oriented_circle_on_sphere_2;

  /// @}

  /// \name Constructions
  ///
  /// @{

  /// Construction object. Must provide the operator:
  ///
  /// `Arc_on_sphere_2 operator()(Point_on_sphere_2 p, Point_on_sphere_2 q)`
  ///
  /// which introduces an arc of great circle, with source `p` and target `q`.
  ///
  /// The circular arc constructed from a source, and a target, is defined as the set of points
  /// of the great circle that lie between the source `p` and the target `q`,
  /// when traversing the great circle counterclockwise seen from the side of the plane
  /// of the great circle pointed by its <I>positive</I> normal vectors.
  ///
  /// In this definition, we say that a normal vector \f$ (a,b,c)\f$ is <I>positive</I> if
  /// \f$ (a,b,c)>(0,0,0)\f$ (i.e.\ \f$ (a>0) || (a==0) \&\& (b>0) || (a==0)\&\&(b==0)\&\&(c>0)\f$).
  ///
  /// \pre `p` and `q` are not diametrically opposed.
  typedef unspecified_type Construct_arc_on_sphere_2;

  /// Construction object. Must provide the operator:
  ///
  /// `Point_on_sphere_2 operator()(Point_3 p)`,
  ///
  /// which expresses the point `p`, a 3D point living on the sphere,
  /// using the (user-defined) point-on-sphere representation.
  typedef unspecified_type Construct_point_on_sphere_2;

  /// Construction object. Must provide the operator:
  ///
  /// `Point_3 operator()(Point_on_sphere_2 p)`,
  ///
  /// which expresses the point `p` on the sphere as a point in the embedding 3D space.
  typedef unspecified_type Construct_point_3;

  /// Construction object. Must provide the operator:
  ///
  /// `Segment_3 operator()(Point_3 p, Point_3 q)`,
  ///
  /// which constructs a segment from two 3D points.
  typedef unspecified_type Construct_segment_3;

  /// Construction object. Must provide the operator:
  ///
  /// `Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r)`,
  ///
  /// which constructs a triangle from three 3D points.
  typedef unspecified_type Construct_triangle_3;

  /// @}

public:
  /// \name Operations
  ///
  /// The following functions provide access to the predicate and constructor objects.
  ///
  /// @{

  // Predicates

  ///
  Compare_on_sphere_2 compare_on_sphere_2_object();

  ///
  Equal_on_sphere_2 equal_on_sphere_2_object();

  ///
  Orientation_on_sphere_2 orientation_on_sphere_2_object();

  ///
  Side_of_oriented_circle_on_sphere_2 side_of_oriented_circle_on_sphere_2_object();

  ///
  Collinear_are_strictly_ordered_on_great_circle_2 collinear_are_strictly_ordered_on_great_circle_2_object();

  // Constructions

  ///
  Construct_arc_on_sphere_2 construct_arc_on_sphere_2_object();

  ///
  Construct_point_on_sphere_2 construct_point_on_sphere_2_object();

  ///
  Construct_point_3 construct_point_3_object();

  ///
  Construct_segment_3 construct_segment_3_object();

  ///
  Construct_triangle_3 construct_triangle_3_object();

  /// @}

public:
  /// \name Configuration of the spherical domain
  ///
  /// @{

  /// sets the center of the sphere.
  ///
  /// \note This function is meant to be used only by the triangulation class as modifying the domain
  ///       requires clearing the triangulation. Users can change the domain using
  ///       `CGAL::Triangulation_on_sphere_2::set_center_and_radius()`.
  void set_center(Point_3 center);

  /// returns the center of the sphere.
  Point_3 center() const;

  /// sets the radius of the sphere.
  ///
  /// \note This function is meant to be used only by the triangulation class as modifying the domain
  ///       requires clearing the triangulation. Users can change the domain using
  ///       `CGAL::Triangulation_on_sphere_2::set_center_and_radius()`.
  void set_radius(FT radius);

  /// returns the radius of the sphere.
  FT radius() const;

  /// @}

public:
  /// \name Precision predicates
  ///
  /// When a kernel does not offer arbitrary precision - or at least an exact representation of
  /// algebraic coordinates, which is usually the case for the sake of computational speed,
  /// then most points that are theoretically on the sphere are not in practice exactly on the sphere.
  /// When this is the case, it cannot be guaranteed that all points are in a convex position,
  /// which can then result in points being hidden upon insertion, see \cgalCite{cgal:ccplr-redtp-10}.
  /// To ensure that no point is hidden, a trick is to enforce a gap between points (\cgalCite{cgal:ccplr-redtp-10},
  /// Lemma 4.1). This gap is based on a fixed maximal allowed distance \f$ \delta \f$ between a point
  /// and the theoretical sphere. The gap must then be at least \f$ 2 \sqrt{R\delta} \f$, where \f$ R \f$
  /// is the radius of the sphere.
  ///
  /// The following two predicates serve to check if the point is on the sphere, that is
  /// if for a model-defined \f$ \delta \f$ the distance to the sphere is smaller
  /// than \f$ \delta \f$ and if the distance between two points is greater than \f$ 2 \sqrt{R\delta} \f$.
  /// It is also of course possible to construct traits classes with arbitrary precision or
  /// with a representation that ensures that points are exactly on the sphere and in this case,
  /// \f$ \delta = 0 \f$ can be chosen and the predicates below are trivial.
  ///
  /// \sa `CGAL::Delaunay_triangulation_on_sphere_traits_2`
  ///
  /// @{

  /// returns whether the point `p` is on the sphere or not.
  ///
  /// \note A point that is not on the sphere will not be inserted in the triangulation.
  bool is_on_sphere(Point_on_sphere_2 p);

  /// returns whether the points `p` and `q` are too close (see
  /// <a href="https://doc.cgal.org/latest/Triangulation_on_sphere_2/classTriangulationOnSphereTraits__2.html#amgrpbc40b59d2ea2be0fc60516e5f7bb0f92">
  /// Precision predicates</a>).
  bool are_points_too_close(Point_on_sphere_2 p, Point_on_sphere_2 q);

  /// @}
};
