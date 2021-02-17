namespace CGAL{

/*!
 \ingroup PkgTriangulation2TriangulationClasses

  The class `Triangulation_2_projection_traits_3` is an extension of the `Projection_traits_xy_3`,
  `Projection_traits_xz_3` and `Projection_traits_yz_3` traits classes, that enables
  using 2D algorithms on the projections of 3D data on any plane.

 \tparam K must be a model of `Kernel`

 \note Internal constructions (projections) are used in the predicate and
construction functors of this class. If `K` is a model of `Kernel` providing exact
constructions or if `K` is a `CGAL::Filtered_kernel` (such as for
`CGAL::Exact_predicates_inexact_constructions_kernel`), this class automatically
provides exact predicates.

\cgalModels `TriangulationTraits_2`
\cgalModels `DelaunayTriangulationTraits_2`
\cgalModels `ConstrainedTriangulationTraits_2`

 */
template <class K>
class Triangulation_2_projection_traits_3
{
public:
  ///\name Types
  ///@{

  //!
  typedef typename K::Point_3     Point_2;

  ///@}


  /// \name Functors
  /// The functors provided by this class are those listed in the
  /// concepts. The functors operate on the 2D projection of their
  /// arguments. They come with preconditions that projections of the
  /// arguments are non-degenerate, eg. a line segment does not project
  /// on a single point, two points do not project on the same point, etc.
  /// The following functor is an addition to the concepts.
  /// @{

  /*!
  \brief projects a 3D point on a plane oriented by the traits.

  A construction object.
  Provides the constructor:

  `Projection_to_plane(const Point_2& plane_point, const Triangulation_2_projection_traits_3& tr)`

  and the operator:

 `Point_2 operator()(const Point_2& point) const`

  which returns the projection of the (3D) point `point` onto the plane passing through
  the point 'plane_point' and with orthogonal vector `normal`, the vector given to the traits constructor.
  */
  typedef unspecified_type Projection_to_plane;
  ///@}

  /// \name Creation
  ///@{

  /*!
   * \brief Constructor
   *
   * \param normal a vector orthogonal to the projection plane.
   */
  Triangulation_2_projection_traits_3(const typename K::Vector_3& normal);

  ///@}
};
} // end namespace CGAL
