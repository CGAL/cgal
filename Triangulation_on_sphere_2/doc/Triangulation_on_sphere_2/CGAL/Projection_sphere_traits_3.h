namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

The class `Projection_sphere_traits_3` is a model of the concept `DelaunayTriangulationOnSphereTraits_2`.

It implements the `Point_on_sphere_2` type as a custom point type by

\tparam K a kernel type; must be a model of `Kernel`

*/
template< typename K >
class Projection_sphere_traits_3
{
public:
  /// The field number type
  typedef typename K::FT                            FT;

  ///
  typedef unspecified_type                          Point_on_sphere_2;

  ///
  typedef typename K::Point_3                       Point_3;

  ///
  typedef typename K::Segment_3                     Segment_3;

  ///
  typedef unspecified_type                          Compare_on_sphere_2;

  ///
  typedef unspecified_type                          Construct_circumcenter_on_sphere_2;

  ///
  typedef typename K::Construct_segment_3           Construct_segment_3;

  /// Two points `p` and `q` are here equal if the rays starting at the center and passing
  /// through `p` and `q` are equal
  typedef unspecified_type                          Equal_on_sphere_2;

  ///
  typedef unspecified_type                          Inside_cone_2;

  ///
  typedef unspecified_type                          Orientation_3;

  ///
  typedef unspecified_type                          Orientation_on_sphere_2;

public:
  /// Due to its representation, any point is exactly on the sphere, and this function
  /// always returns `true`.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// Since there is no need to ensure separations of the points because representation
  /// of the points is exact (see also \cgalCite{cgal:ccplr-redtp-10}),
  /// this function simply returns `false` for any input.
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;
};

} // namespace CGAL
