
/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The parameterized plane shape that is fitted.

\cgalRefines `Proxy`
\cgalHasModel `PlaneProxy`

*/
class PlaneProxy {
public:
  /*! 3D point type
   * It must be default constructible, and can be constructed from 3 objects of type `FT`.
   * `bool operator<(Point_3, Point_3)` to lexicographically compare two points must be available.
   * Access to Cartesian coordinates must be possible using `Point_3::x()`, `Point_3::y(), Point_3::z()` and
   * `FT operator[](int i)` with  `0 <= i < 3`.
   *
   * There must be a specialization of `CGAL::Kernel_traits` such that
   * `CGAL::Kernel_traits<Point_3>::%Kernel` is a model implementing this concept.
   */
  typedef unspecified_type Point_3;

  /// 3D vector type
  typedef unspecified_type Vector_3;
  /// 3D plane type
  typedef unspecified_type Plane_3;

  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;

  /// @name Data members
  /// @{

  /// data member to describe the plane proxy center position.
  Point_3 center;
  /// data member to describe the plane proxy normal.
  Vector_3 normal;
  /// data member to describe the proxy seed.
  facet_descriptor seed;
  /// data member to describe the fitted plane.
  Plane_3 fit_plane;

  /// }
};

