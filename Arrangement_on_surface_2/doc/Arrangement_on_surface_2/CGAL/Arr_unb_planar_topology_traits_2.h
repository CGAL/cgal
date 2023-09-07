namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Ref
 *
 * \anchor arr_ref_unb_planar_topology_traits
 *
 * This class handles the topology for arrangements of unb curves
 * embedded in the plane.
 *
 * The `Arr_unb_planar_topology_traits_2` template has two parameters:
 * <UL>
 * <LI>The `GeometryTraits_2` template-parameter should be instantiated with
 * a model of the `ArrangementBasicTraits_2` concept. The traits
 * class defines the types of \f$x\f$-monotone curves and two-dimensional
 * points, namely `ArrangementBasicTraits_2::X_monotone_curve_2` and
 * `ArrangementBasicTraits_2::Point_2`,
 *   respectively, and supports basic geometric predicates on them.
 * <LI>The `Dcel` template-parameter should be instantiated with
 * a class that is a model of the `ArrangementDcel` concept. The
 * value of this parameter is by default
 * `Arr_default_dcel<Traits>`.
 * </UL>
 *
 * \cgalModels{ArrangementBasicTopologyTraits}
 *
 * \sa `Arr_default_dcel<Traits>`
 * \sa `CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel,x,y>`
 */
template <typename GeometryTraits_2,
          typename Dcel = Arr_default_dcel<GeometryTraits_2> >
class Arr_unb_planar_topology_traits_2 {
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

  /*! Default constructor. */
  Arr_unb_planar_topology_traits_2();

  /*! Constructor from a geometry-traits object.
   * \param traits the traits.
   */
  Arr_unb_planar_topology_traits_2(const GeometryTraits_2* traits);

  /// @}

  /// \name Accessors
  /// @{

  /*! Obtain the DCEL (const version). */
  const Dcel& dcel() const;

  /*! Obtain the DCEL (non-const version). */
  Dcel& dcel();

  /*! Obtain the unbounded face (const version). */
  const Face* unbounded_face() const;

  /*! Obtain the unbounded face (non-const version). */
  Face* unbounded_face();

  /// @}

  /// \name Modifiers
  /// @{
  /// @}
};

}
