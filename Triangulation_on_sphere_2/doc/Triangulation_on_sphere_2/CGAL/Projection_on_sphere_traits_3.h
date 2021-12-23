namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Projection_on_sphere_traits_3` is a model of the concept `DelaunayTriangulationOnSphereTraits_2`.

It implements the `Point_on_sphere_2` type as a custom point type which represents the projection
of a point living in the 3D Euclidean space onto the sphere along the segment between said point
and the center of the sphere.

\tparam LK a linear kernel type; must be a model of `Kernel`.
\tparam SK a spherical kernel type; must be a model of `SphericalKernel`.

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

\sa `CGAL::Delaunay_triangulation_on_sphere_traits_2`
*/
template <typename LK,
          typename SK = CGAL::Spherical_kernel_3<
                          LK, CGAL::Algebraic_kernel_for_spheres_2_3<typename LK::FT> > >
class Projection_on_sphere_traits_3
{
public:
  /// The field number type
  typedef typename LK::FT                           FT;

  /// The point on the sphere type
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

  /// Points are equal if they have the same projection onto the sphere.
  typedef unspecified_type                          Equal_on_sphere_2;

  /// @}

public:
  /// \name Precision predicates
  ///
  /// @{

  /// Due to the chosen point representation, any point is theoretically on the sphere, and this function
  /// always returns `true`.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// returns `false` if `LK` can represent algebraic coordinates, or whether the distance
  /// between `p` and `q` is lower than \f$ 2 \sqrt{R\delta} \f$ otherwise (see the traits class
  /// `CGAL::Delaunay_triangulation_on_sphere_traits_2`).
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;

  /// @}
};

} // namespace CGAL
