namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Delaunay_triangulation_on_sphere_traits_2` is a model
of the concept `DelaunayTriangulationOnSphereTraits_2`.

The `Point_on_sphere_2` type is implemented as a kernel `Point_3` type.

If the kernel template parameter `LK` does not enable exact representation of points on sphere
(i.e., at least a mean to represent algebraic coordinates), then it cannot be guaranteed
that all points on the sphere are in a convex position and thus some points might be hidden upon insertion.
It is possible to ensure that no point can be hidden by enforcing a tiny gap between points \cgalCite{cgal:ccplr-redtp-10}.
In Lemma 4.1 of this publication, it is in particular proven that if points lie within a distance
\f$ \delta \f$ of the sphere, no point is hidden as long as points are separated by at least \f$ 2 \sqrt{R\delta} \f$
with \f$ R \f$ the radius of the sphere.

Thus, if `LK` offers exact representation, then \f$ \delta \f$ is set to \f$ 0 \f$.
Otherwise, \f$ \delta \f$ is set to the maximal distance between two consecutive `LK::FT`
for the coordinates of a point that is (theoretically) on the sphere.
This bound \f$ \delta \f$ is then used in the functions `is_on_sphere()` and `are_points_too_close()`
using the relation above to ensure that a point being inserted is either marked as "too close"
(see \link CGAL::Triangulation_on_sphere_2::Locate_type `CGAL::Triangulation_on_sphere_2<Traits,TDS>::Locate_type` \endlink)
and thus not inserted, or guaranteed to not be hidden upon insertion.

\tparam LK a linear kernel type; it must be a model of `Kernel`.
\tparam SK a spherical kernel type; it must be a model of `SphericalKernel` refining `LK`.

\cgalModels{DelaunayTriangulationOnSphereTraits_2}

\sa `CGAL::Projection_on_sphere_traits_3`
*/
template <typename LK,
          typename SK = CGAL::Spherical_kernel_3<
                          LK, CGAL::Algebraic_kernel_for_spheres_2_3<typename LK::FT> > >
class Delaunay_triangulation_on_sphere_traits_2
{
public:
  /// The field number type
  typedef typename LK::FT                           FT;

  ///
  typedef typename LK::Point_3                      Point_on_sphere_2;

  /// An arc of a great circle, used to represent a curved segment on the sphere (Voronoi or Delaunay edge).
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

  /// Internally uses `LK::Coplanar_orientation_3`
  typedef unspecified_type                          Collinear_are_strictly_ordered_on_great_circle_2;

  ///
  typedef typename LK::Compare_xyz_3                Compare_on_sphere_2;

  /// If the kernel cannot represent algebraic coordinates exactly, there is a tolerance
  /// around the sphere, and thus different points of \f$ \mathbb{R}^3\f$ can actually
  /// correspond to the same point on the sphere.
  /// This functor checks if two points project onto the same point on the sphere.
  typedef unspecified_type                          Equal_on_sphere_2;

  ///
  typedef typename LK::Orientation_3                Side_of_oriented_circle_on_sphere_2;

  /// Internally uses `LK::Orientation_3`
  typedef unspecified_type                          Orientation_on_sphere_2;

  /// @}

  /// \name Constructions
  ///
  /// @{

  /// Internally uses `SK::Construct_circular_arc_3`
  typedef typename unspecified_type                 Construct_arc_on_sphere_2;

  ///
  typedef typename LK::Construct_circumcenter_3     Construct_circumcenter_3;

  ///
  typedef typename LK::Construct_point_3            Construct_point_3;

  ///
  typedef typename LK::Construct_segment_3          Construct_segment_3;

  ///
  typedef typename LK::Construct_triangle_3         Construct_triangle_3;

  /// @}

public:
  /// \name Precision predicates
  ///
  /// @{

  /// returns whether `p` is exactly on the sphere if `LK` can represent algebraic coordinates,
  /// or whether `p` is within an automatically computed small distance otherwise.
  bool is_on_sphere(const Point_on_sphere_2& p) const;

  /// returns `false` if `LK` can represent algebraic coordinates, or whether the distance
  /// between `p` and `q` is lower than \f$ 2 \sqrt{R\delta} \f$ otherwise.
  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const;

  /// @}
};

} // namespace CGAL
