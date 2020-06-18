namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

\cgalModels `DelaunayTriangulationOnSphereTraits_2`

The class `Delaunay_triangulation_sphere_traits_2` is a model
of the concept `DelaunayTriangulationOnSphereTraits_2`. It implements the `Point_on_sphere_2` type
as a kernel's `Point_3` type.

If the kernel template parameter `K` does not enable exact representation of points on sphere
(i.e. at least a mean to represent algebraic coordinates), then it cannot be guaranteed
that all points on the sphere are in a convex position and thus some points might be hidden upon insertion.
It is possible to ensure that no point can be hidden by enforcing a tiny gap between points \cgalCite{cgal:ccplr-redtp-10}.
In Lemma 4.1 of this publication, it is in particular proven that if points lie within a distance
\f$ \delta \f$ of the sphere, then as long as points are separated by at least \f$ 2 \sqrt{R\delta} \f$
with \f$ R \f$ the radius of the sphere, then no point is hidden.

Thus, if `K` offers exact representation, then \f$ \delta \f$ is set to \f$ 0 \f$.
Otherwise, \f$ \delta \f$ is set to the maximal distance between two consecutive `K::FT`
for the coordinates of a point that is (theoretically) on the sphere.
This bound \f$ \delta \f$ is then used in the functions `is_on_sphere()` and `are_points_too_close()`
using the relation above to ensure that a point being inserted is either marked as "too close"
(see \link CGAL::Triangulation_on_sphere_2::Locate_type `CGAL::Triangulation_on_sphere_2<Traits,TDS>::Locate_type` \endlink)
and not inserted or guaranteed to not be hidden upon insertion.

\tparam K a kernel type; must be a model of `Kernel`
\tparam SK a spherical kernel type; must be a model of `SphericalKernel`

\sa `CGAL::Geographical_coordinates_traits_2`
\sa `CGAL::Projection_sphere_traits_3`
*/
template <typename K,
          typename SK = CGAL::Spherical_kernel_3<
                          K, CGAL::Algebraic_kernel_for_spheres_2_3<typename K::FT> > >
class Delaunay_triangulation_sphere_traits_2
{
public:
  /// The field number type
  typedef typename K::FT                            FT;

  ///
  typedef typename K::Point_3                       Point_on_sphere_2;

  /// An arc of a great circle, used to represent a curved segment on the sphere (Voronoi or Delaunay edge).
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

  ///
  typedef typename K::Compare_xyz_3                 Compare_on_sphere_2;

  /// If the kernel cannot represent algeabric coordinates exactly, there is a tolerance
  /// around the sphere, and thus different points can actually be the same point.
  /// This particular equality functor checks if both query points are on the sphere and
  /// are aligned (and on the same side) with the center of the sphere.
  typedef unspecified_type                          Equal_on_sphere_2;

  /// Internally uses a `K::Coplanar_orientation_3`
  typedef unspecified_type                          Collinear_are_strictly_ordered_on_great_circle_2;

  ///
  typedef typename K::Orientation_3                 Side_of_oriented_circle_2;

  /// Internally uses `Orientation_3`
  typedef unspecified_type                          Orientation_on_sphere_2;

  /// @}

  /// \name Constructions
  ///
  /// @{

  /// Internally uses `SK::Construct_circular_arc_3`
  typedef typename unspecified_type                 Construct_arc_on_sphere_2;

  ///
  typedef typename unspecified_type                 Construct_circumcenter_on_sphere_2;

  ///
  typedef typename K::Construct_point_3             Construct_point_3;

  ///
  typedef typename K::Construct_segment_3           Construct_segment_3;

  ///
  typedef typename K::Construct_triangle_3          Construct_triangle_3;

  /// @}

public:
  /// \name Precision predicates
  ///
  /// @{

  /// Returns whether `p` is exactly on the sphere if `K` can represent algebraic coordinates,
  /// or whether `p` is within an automatically computed small distance otherwise.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// Returns `false` if `K` can represent algeabric coordinates, or whether the distance
  /// between `p` and `q` is lower than \f$ 2 \sqrt{R\delta} \f$ otherwise.
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;

  /// @}
};

} // namespace CGAL
