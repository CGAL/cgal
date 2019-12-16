namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

The class `Delaunay_triangulation_sphere_traits_2` is a model
of the concept `DelaunayTriangulationOnSphereTraits_2`. It implements the `Point_on_sphere` type
as a kernel's `Point_3` type.

If the kernel template parameter `K` does not enable exact representation of points on sphere
(i.e. at least a mean to represent algebraic coordinates), then it cannot be guaranteed
that all points on the sphere are in a convex position and thus some points might be hidden upon insertion.
It is possible to ensure that no point can be hidden by enforcing a tiny gap between points \cgalCite{cgal:ccplr-redtp-10}.
In Lemma 4.1 of this publication, it is in particular proven that if points lie within a distance
\f$ \delta \f$ of the sphere, then as long as points are separated by at least \f$ 2 \sqrt{R\delta} \f$
with \f$ R \f$ the radius of the sphere, then no point is hidden.

Thus, if `K` offers exact representation, then \f$ \delta \f$ is set to \f$ 0 \f$.
Otherwise,  \f$ \delta \f$ is set to the maximal distance between two consecutive `K::FT`
for the coordinates of a point that is (theoretically) on the sphere.
This bound  \f$ \delta \f$ is then used in the functions `is_on_sphere()` and `are_points_too_close()`
using the relation above to ensure that a point being inserted is either  marked as "too close"
(see \link CGAL::Triangulation_on_sphere_2::Locate_type `CGAL::Triangulation_on_sphere_2<Traits,TDS>::Locate_type` \endlink)
and not inserted or guaranteed to not be hidden upon insertion.

@todo Right now it takes 2^-25*radius, but we could (should) get the smallest distance between
      two doubles by using std::nextafter at 'center.x() + radius'. Then we can deduce whatever is needed
      for kernels that have a simple type. However, what to do for kernels that have weird types,
      e.g. someone that wants to run this with simple_cartesian<custom_FT>? How to compute the precision?
      Just leave it to the user to implement his own class then by throwing an error "couldn't
      deduce necessary precision"?

\tparam K a kernel type; must be a model of `Kernel`

\sa `CGAL::Geographical_coordinates_traits_2`
\sa `CGAL::Projection_sphere_traits_3`
*/
template< typename K >
class Delaunay_triangulation_sphere_traits_2
{
public:
  /// The field number type
  typedef typename K::FT                            FT;

  /// For this model, the point on the sphere is simply a 3D point
  typedef typename K::Point_3                       Point_on_sphere;

  ///
  typedef typename K::Point_3                       Point_3;

  ///
  typedef typename K::Segment_3                     Segment_3;

  ///
  typedef typename K::Compare_xyz_3                 Compare_on_sphere_2;

  ///
  typedef typename K::Construct_circumcenter_3      Construct_circumcenter_3;

  ///
  typedef unspecified_type                          Construct_point_3;

  ///
  typedef typename K::Construct_segment_3           Construct_segment_3;

  /// Internally uses `K::Collinear_3`
  typedef unspecified_type                          Equal_on_sphere_2;

  /// Internally uses a `K::Coplanar_orientation_3`
  typedef unspecified_type                          Inside_cone_2;

  ///
  typedef typename K::Orientation_3                 Orientation_3;

  /// Internally uses `Orientation_3`
  typedef unspecified_type                          Orientation_on_sphere_2;

public:
  /// Returns whether `p` is exactly on the sphere if `K` can represent algebraic coordinates,
  /// or whether `p` is within an automatically computed small distance otherwise.
  bool is_on_sphere(const Point_on_sphere& p) const;

  /// Returns `false` if `K` can represent algeabric coordinates, and whether the distance
  /// between `p` and `q` is lower than \f$ 2 \sqrt{R\delta} \f$ otherwise.
  bool are_points_too_close(const Point_on_sphere& p, const Point_on_sphere& q) const;
};

} // namespace CGAL
