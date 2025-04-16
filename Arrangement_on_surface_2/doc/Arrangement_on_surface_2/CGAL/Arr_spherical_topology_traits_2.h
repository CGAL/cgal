namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Ref
 *
 * \anchor arr_ref_spherical_topology_traits
 *
 * This class handles the topology for arrangements of great spherical
 * arcs on the sphere embedded on 2D parametric surdace.
 *
 * The `Arr_spherical_topology_traits_2` template has two parameters:
 * <UL>
 * <LI>The `GeometryTraits_2` template-parameter should be substituted by
 * a model of the `AosBasicTraits_2` concept. The traits
 * class defines the types of \f$x\f$-monotone curves and two-dimensional
 * points, namely `AosBasicTraits_2::X_monotone_curve_2` and
 * `AosBasicTraits_2::Point_2`,
 *   respectively, and supports basic geometric predicates on them.
 * <LI>The `Dcel` template-parameter should be substituted by
 * a class that is a model of the `AosDcel` concept. The
 * value of this parameter is by default
 * `Arr_default_dcel<Traits>`.
 * </UL>
 *
 * \cgalModels{AosBasicTopologyTraits}
 *
 * \sa `Arr_default_dcel<Traits>`
 * \sa `CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel,x,y>`
 */
template <typename GeometryTraits_2,
          typename Dcel = Arr_default_dcel<GeometryTraits_2> >
class Arr_spherical_topology_traits_2 {
public:
  /// \name Types
  /// @{

  typedef typename GeometryTraits_2::Point_2            Point_2;
  typedef typename GeometryTraits_2::X_monotone_curve_2 X_monotone_curve_2;

  typedef typename Dcel::Size                           Size;
  typedef typename Dcel::Vertex                         Vertex;
  typedef typename Dcel::Halfedge                       Halfedge;
  typedef typename Dcel::Face                           Face;
  typedef typename Dcel::Outer_ccb                      Outer_ccb;
  typedef typename Dcel::Inner_ccb                      Inner_ccb;
  typedef typename Dcel::Isolated_vertex                Isolated_vertex;

  /// @}

  /// \name Creation
  /// @{

  /*! constructs default */
  Arr_spherical_topology_traits_2();

  /*! constructs from a geometry-traits object.
   * \param traits the traits.
   */
  Arr_spherical_topology_traits_2(const GeometryTraits_2* traits);

  /// @}

  /// \name Accessors
  /// @{

  /*! obtains the \dcel (const version). */
  const Dcel& dcel() const;

  /*! obtains the \dcel (non-const version). */
  Dcel& dcel();

  /*! obtains the spherical face (const version). */
  const Face* spherical_face() const;

  /*! obtains the spherical face (non-const version). */
  Face* spherical_face();

  /*! obtains the south pole (const version). */
  const Vertex* south_pole() const;

  /*! obtains the south pole (non-const version). */
  Vertex* south_pole();

  /*! obtains the north pole (const version). */
  const Vertex* north_pole() const;

  /*! obtains the north pole (non-const version). */
  Vertex* north_pole();

  /*! obtains a vertex on the line of discontinuity that corresponds to
   *  the given point (or return NULL if no such vertex exists).
   */
  const Vertex* discontinuity_vertex(const Point_2& pt) const;

  /*! obtains a vertex on the line of discontinuity that corresponds to
   *  the given point (or return NULL if no such vertex exists).
   */
  Vertex* discontinuity_vertex(const Point_2& pt);

  /// @}

  /// \name Modifiers
  /// @{
  /// @}
};

}
