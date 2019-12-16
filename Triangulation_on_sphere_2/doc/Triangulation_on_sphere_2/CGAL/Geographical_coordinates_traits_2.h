namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

The class `Geographical_coordinates_traits_2` is a model
of the concept `DelaunayTriangulationOnSphereTraits_2`. It implements the `Point_on_sphere` type
using a custom coordinate type made.

\tparam K a kernel type; must be a model of `Kernel`

\sa `CGAL::Delaunay_triangulation_sphere_traits_2`
\sa `CGAL::Projection_sphere_traits_3`
*/
template< typename K >
class Geographical_coordinates_traits_2
{
public:
  /// The field number type.
  typedef typename K::FT                            FT;

  /// A pair of lattitude and longitude values.
  typedef Geographical_coordinates<K>               Point_on_sphere;

  ///
  typedef typename K::Point_3                       Point_3;

  ///
  typedef typename K::Segment_3                     Segment_3;

  /// Uses a 2-dimensional lexicographical order to create a strict total order on the sphere.
  typedef unspecified_type                          Compare_on_sphere_2;

  ///
  typedef unspecified_type                          Construct_circumcenter_3;

  ///
  typedef unspecified_type                          Construct_point_3;

  ///
  typedef typename K::Construct_segment_3           Construct_segment_3;

  /// Two points are equal if their two coordinates (lattitude and longitude) are equal.
  typedef unspecified_type                          Equal_on_sphere_2;

  ///
  typedef unspecified_type                          Inside_cone_2;

  ///
  typedef unspecified_type                          Orientation_3;

  ///
  typedef unspecified_type                          Orientation_on_sphere_2;

public:
  /// Due to their representation, points are always exactly on the sphere, and consequently
  /// this function simply returns `true` for any input.
  bool is_on_sphere(const Point_on_sphere& p) const;

  /// Since there is no need to ensure separations of the points because representation
  /// of the points is exact (see also \cgalCite{cgal:ccplr-redtp-10}),
  /// this function simply returns `false` for any input.
  bool are_points_too_close(const Point_on_sphere& p, const Point_on_sphere& q) const;
};

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

This class represents coordinates of the Geographical Coordinates System,
that is a pair of two values representing a longitude and a lattitude.

@todo should this be in radians?

\tparam K a kernel type; must be a model of `Kernel`

*/
template< typename K >
class Geographical_coordinates
{
public:
  /// The field number type
  typedef typename K::FT                            FT;

  ///
  typedef FT                                        Lattitude;

  ///
  typedef FT                                        Longitude;

  /// %Default constructor. Creates a point at coordinates `(0, 0)`.
  Geographical_coordinates();

  /// Construct a point on the sphere at coordinates `(la, lo)`.
  ///
  /// \pre `la` is within `[-90; 90[` and `lo` is within `[-180; 180[`.
  Geographical_coordinates(const Lattitude la, const Longitude lo);
};

} // namespace CGAL
