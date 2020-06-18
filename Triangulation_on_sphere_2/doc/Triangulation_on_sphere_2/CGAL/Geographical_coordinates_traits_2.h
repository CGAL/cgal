namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

The class `Geographical_coordinates_traits_2` is a model
of the concept `DelaunayTriangulationOnSphereTraits_2`. It implements the `Point_on_sphere_2` type
as a pair of coordinates representing the latitude and the longitude of the point on the sphere.

\tparam K a kernel type; must be a model of `Kernel`
\tparam SK a spherical kernel type; must be a model of `SphericalKernel`

\sa `CGAL::Delaunay_triangulation_sphere_traits_2`
\sa `CGAL::Projection_sphere_traits_3`
*/
template <typename K,
          typename SK = CGAL::Spherical_kernel_3<
                          K, CGAL::Algebraic_kernel_for_spheres_2_3<typename K::FT> > >
class Geographical_coordinates_traits_2
{
public:
  /// The field number type.
  typedef typename K::FT                            FT;

  /// A pair of latitude and longitude values.
  typedef Geographical_coordinates<K>               Point_on_sphere_2;

  ///
  typedef typename SK::Circular_arc_3               Arc_on_sphere_2;

  ///
  typedef typename K::Point_3                       Point_3;

  ///
  typedef typename K::Segment_3                     Segment_3;

  ///
  typedef typename K::Triangle_3                    Triangle_3;

  /// \name Predicates
  ///
  /// @{

  /// The 2-dimensional lexicographical order is used to create a strict total order on the sphere.
  typedef unspecified_type                          Compare_on_sphere_2;

  /// Two points are equal if their two coordinates (latitude and longitude) are equal.
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
  typedef typename unspecified_type                 Construct_arc_on_sphere_2;

  ///
  typedef typename unspecified_type                 Construct_circumcenter_on_sphere_2;

  /// Converts points from the latitude/longitude system to the 3D Euclidean system
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

  /// Due to their representation, points are always exactly on the sphere, and consequently
  /// this function simply returns `true` for any input.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// Since there is no need to ensure separations of the points because representation
  /// of the points is exact (see also \cgalCite{cgal:ccplr-redtp-10}),
  /// this function simply returns `false` for any input.
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;

  /// @}
};

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

This class represents coordinates of the Geographical Coordinates System,
that is a pair of two values representing a longitude and a latitude.

\tparam K a kernel type; must be a model of `Kernel`

*/
template< typename K >
class Geographical_coordinates
{
public:
  /// The field number type
  typedef typename K::FT                            FT;

  ///
  typedef FT                                        latitude;

  ///
  typedef FT                                        Longitude;

  /// %Default constructor. Creates a point at coordinates `(0, 0)`.
  Geographical_coordinates();

  /// Construct a point on the sphere at coordinates `(la, lo)`.
  ///
  /// \pre `la` is within `[-90; 90[` and `lo` is within `[-180; 180[`.
  Geographical_coordinates(const latitude la, const Longitude lo);
};

} // namespace CGAL
