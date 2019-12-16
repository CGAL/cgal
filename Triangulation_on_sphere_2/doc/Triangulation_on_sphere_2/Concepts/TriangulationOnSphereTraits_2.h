/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

The concept `TriangulationOnSphereTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Triangulation_on_sphere_2<Traits, Tds>`.
This concept provides the types of the geometric primitives used in the
triangulation and some function object types for the required predicates on those primitives.

In particular, the traits class is expected to contain information about the sphere (center and radius)
on which the points of the triangulation lie, as well as to provide means to change this information.

\cgalHasModel `CGAL::Delaunay_triangulation_sphere_traits_2<K>`
\cgalHasModel `CGAL::Geographical_coordinates_traits_2<K>`
\cgalHasModel `CGAL::Projection_sphere_traits_3<K>`

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
  typedef unspecified_type Point_on_sphere;

  /*!
  The 3D point type.
  */
  typedef unspecified_type Point_3;

public:
  // Predicates and constructions

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere p, Point_on_sphere q)`
  ///
  /// which returns `true` or `false` based on some (arbitrary) strict total order of points.
  typedef unspecified_type Compare_on_sphere_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere p, Point_on_sphere q)`
  ///
  /// which returns `true` if and only if `p` and `q` are equal.
  typedef unspecified_type Equal_on_sphere_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `bool operator()(Point_on_sphere p, Point_on_sphere q, Point_on_sphere r)`
  ///
  /// which returns `true` if `r` strictly lies inside the cone defined by the center of the sphere,
  /// `p` and `q`.
  typedef unspecified_type Inside_cone_2;

  /// Predicate object type. Must provide the operator:
  ///
  /// `Orientation operator()(Point_on_sphere p, Point_on_sphere q, Point_on_sphere r, Point_on_sphere s)`
  ///
  /// which returns `CGAL::POSITIVE`, if `s` lies on the positive side of the oriented plane `h`
  /// defined by `p`, `q`, and `r`; returns `CGAL::NEGATIVE` if `s` lies on the negative side of `h`;
  /// and returns `CGAL::COPLANAR` if `s` lies on `h`.
  typedef unspecified_type Orientation_3;

  /// Predicate object type. Must provide the operators:
  ///
  /// `Orientation operator()(Point_on_sphere p, Point_on_sphere q, Point_on_sphere r)`
  ///
  /// which returns `CGAL::POSITIVE`, if `r` lies on the positive side of the oriented plane `h`
  /// defined by the center of the sphere, `p`, and `q`; returns `CGAL::NEGATIVE` if `r`
  /// lies on the negative side of `h`; and returns `CGAL::COPLANAR` if `r` lies on `h`.
  ///
  /// \note This is distinguished from `Orientation_3` because the center of the sphere has type
  /// `Point_3`, which might not be equivalent or convertible to `Point_on_sphere`.
  typedef unspecified_type Orientation_on_sphere_2;

public:
  /// \name Operations
  ///
  /// The following functions give access to the predicate and constructor objects.
  ///
  /// @{

  ///
  Compare_on_sphere_2 compare_on_sphere_2_object();

  ///
  Equal_on_sphere_2 equal_on_sphere_2_object();

  ///
  Inside_cone_2 inside_cone_2_object();

  ///
  Orientation_3 orientation_3_object();

  ///
  Orientation_on_sphere_2 orientation_on_sphere_2_object();

  /// @}

public:
  /// \name Configuration of the spherical domain
  ///
  /// @{

  /// Sets the center of the sphere.
  void set_center(Point_3);

  /// Returns the center of the sphere.
  Point_3 center();

  /// Sets the radius of the sphere.
  void set_radius(FT radius);

  /// Returns the radius of the sphere.
  FT radius();

  /// @}

public:
  /// \name Precision predicates
  ///
  /// When a kernel does not offer arbitrary precision - or at least an exact representation of
  /// algebraic coordinates, which is usually the case for the sake of computational speed,
  /// then most points that are theoretically on the sphere are not in practice exactly in the sphere.
  /// When this is the case, it cannot be guaranteed that all points are in a convex position,
  /// which can then result in points being hidden upon insertion, see \cgalCite{cgal:ccplr-redtp-10}.
  /// To ensure that no point is hidden, a trick is to enforce a gap between points (\cgalCite{cgal:ccplr-redtp-10},
  /// Lemma 4.1). This gap is based on a fixed maximal allowed distance \f$ \delta \f$ between a point
  /// and the theoretical sphere. The gap must then be at least \f$ 2 \sqrt{R\delta} \f$, where \f$ R \f$
  /// is the radius of the sphere.
  ///
  /// The following two predicates serve to check if the point is on the sphere, that is, if for
  /// a model-defined \f$ \delta \f$ the distance to the sphere is smaller than \f$ \delta \f$,
  /// and if the distance between two points is greater than \f$ 2 \sqrt{R\delta} \f$.
  /// It is also of course possible to construct traits classes with arbitrary precision or
  /// with a representation that ensures that points are exactly on the sphere (see for
  /// example `CGAL::Geographical_coordinates_traits_2`) and in this case,
  /// \f$ \delta = 0 \f$ can be chosen and the predicates below are trivial.
  ///
  /// \sa `CGAL::Delaunay_triangulation_sphere_traits_2`
  ///
  /// @{

  /// Returns whether the point `p` is on the sphere or not.
  ///
  /// \note A point that is not on the sphere will not be inserted in the triangulation.
  bool is_on_sphere(Point_on_sphere p);

  /// Returns whether `p` and `q` are too close.
  bool are_points_too_close(Point_on_sphere p, Point_on_sphere q);
};
