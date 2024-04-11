// Undocumented for now
#ifndef DOXYGEN_RUNNING

/*!
\ingroup PkgAlphaWrap3Concepts
\cgalConcept

The concept `AlphaWrapOracle` defines the requirements for an Alpha Wrap <em>Oracle</em>, that is a class
that answers a number of queries over the input of the algorithm.
The oracle is the template parameter of the class `CGAL::Alpha_wraps_3_::Alpha_wrap_3`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Alpha_wraps_3_::Point_set_oracle}
\cgalHasModels{CGAL::Alpha_wraps_3_::Segment_soup_oracle}
\cgalHasModels{CGAL::Alpha_wraps_3_::Triangle_mesh_oracle}
\cgalHasModels{CGAL::Alpha_wraps_3_::Triangle_soup_oracle}
\cgalHasModelsEnd

*/
template <typename GeomTraits>
class AlphaWrapOracle
{
public:
  /// The geometric traits, must be a model of `AlphaWrapOracleTraits_3`
  typedef unspecified_type Geom_traits;

  /// Field type
  typedef Geom_traits::FT FT;

  /// Point type
  typedef Geom_traits::Point_3 Point_3;

  /// Sphere type
  typedef Geom_traits::Ball_3 Ball_3;

  /// Returns the geometric traits
  Geom_traits geom_traits();

  /// Returns an axis-aligned box enclosing the input data.
  CGAL::Bbox_3 bbox();

  /// Returns whether the ball `b` intersects the input data.
  bool do_intersect(Ball_3 b);

  /// Returns whether the tetrahedron `tet` intersects the input data.
  template <typename K>
  bool do_intersect(Tetrahedron_with_outside_info<K> tet);

  /// Returns the intersection `o` of the segment `[p;q]` with the offset isolevel at distance `os`
  /// from the input. In case of multiple intersections, the intersection closest to `p` is returned.
  /// Returns `true` if there is an intersection, and `false` otherwise.
  bool first_intersection(Point_3 p, Point_3 q, Point_3& o, double os);

  /// Returns the point on the input data closest to `q`.
  Point_3 closest_point(Point_3 q);

  /// Returns the smallest squared distance between `q` and the input data.
  FT squared_distance(Point_3 q);
};

#endif // DOXYGEN_RUNNING
