namespace CGAL{

/*!
\ingroup kernel_classes

The class `Projection_traits_3` works similarly to the `Projection_traits_xy_3`,
`Projection_traits_xz_3`, and `Projection_traits_yz_3` traits classes, enabling
the use of 2D algorithms on the projections of 3D data onto an arbitrary plane.

\tparam K must be a model of `Kernel`

\note Internal constructions (projections) are used in the predicate and
construction functors of this class. If `K` is a model of `Kernel` providing exact
constructions or if `K` is a `CGAL::Filtered_kernel` (such as for
`CGAL::Exact_predicates_inexact_constructions_kernel`), this class automatically
provides exact predicates.

\cgalModels{TriangulationTraits_2,DelaunayTriangulationTraits_2,ConstrainedTriangulationTraits_2,PolygonTraits_2,ConformingDelaunayTriangulationTraits_2,Barycentric_coordinates::BarycentricTraits_2}

\sa `CGAL::Projection_traits_xy_3`
\sa `CGAL::Projection_traits_xz_3`
\sa `CGAL::Projection_traits_yz_3`
*/
template <class K>
class Projection_traits_3
{
public:
  /// \name Functors
  /// The functors provided by this class are those listed in the
  /// concepts. The functors operate on the 2D projections of their
  /// arguments. They come with preconditions that projections of the
  /// arguments are non-degenerate, e.g. a line segment does not project
  /// on a single point, two points do not project onto the same point, etc.

  /// \name Creation
  ///@{

  /*!
   * \brief Constructor
   *
   * \param normal a vector orthogonal to the projection plane.
   */
  Projection_traits_3(const typename K::Vector_3& normal);

  ///@}
};

} // end namespace CGAL
