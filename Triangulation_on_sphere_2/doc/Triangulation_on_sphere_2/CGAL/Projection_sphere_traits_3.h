namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Projection_sphere_traits_3` is a model of the concept `DelaunayTriangulationOnSphereTraits_2`.

It implements the `Point_on_sphere_2` type as a custom point type which represents the projection
of a point living in the 3D Euclidean space onto the sphere along the segment between said point
and the center of the sphere.

\tparam K a linear kernel type; must be a model of `Kernel`
\tparam SK a spherical kernel type; must be a model of `SphericalKernel`

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

\sa `CGAL::Delaunay_triangulation_sphere_traits_2`
\sa `CGAL::Geographical_coordinates_traits_2`
*/
template <typename LK,
          typename SK = CGAL::Spherical_kernel_3<
                          K, CGAL::Algebraic_kernel_for_spheres_2_3<typename LK::FT> > >
class Projection_sphere_traits_3
{
public:
  /// The field number type.
  typedef typename LK::FT                           FT;

  /// An internal point type representing the projected point.
  typedef unspecified_type                          Point_on_sphere_2;

  ///
  typedef typename SK::Circular_arc_3               Arc_on_sphere_2;

  ///
  typedef typename LK::Point_3                      Point_3;

  ///
  typedef typename LK::Segment_3                    Segment_3;

  ///
  typedef typename LK::Triangle_3                   Triangle_3;

  /// \name Predicates
  ///
  /// @{

  ///
  typedef unspecified_type                          Compare_on_sphere_2;

  /// Points are equal if they have the same projection onto the sphere
  typedef unspecified_type                          Equal_on_sphere_2;

  ///
  typedef unspecified_type                          Collinear_are_strictly_ordered_on_great_circle_2;

  ///
  typedef unspecified_type                          Side_of_oriented_circle_2;

  ///
  typedef unspecified_type                          Orientation_on_sphere_2;

  /// @}

  /// \name Constructions
  ///
  /// @{

  ///
  typedef unspecified_type                          Compare_on_sphere_2;

  ///
  typedef unspecified_type                          Construct_circumcenter_on_sphere_2;

  ///
  typedef unspecified_type                          Construct_point_3;

  ///
  typedef unspecified_type                          Construct_segment_3;

  ///
  typedef unspecified_type                          Construct_triangle_3;

  /// @}

public:
  /// \name Precision predicates
  ///
  /// @{

  /// Due to its representation, any point is exactly on the sphere, and this function
  /// always returns `true`.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// Since there is no need to ensure separations of the points because the representation
  /// of the points is exact (see also \cgalCite{cgal:ccplr-redtp-10}),
  /// this function simply returns `false` for any input.
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;

  /// @}
};

} // namespace CGAL
