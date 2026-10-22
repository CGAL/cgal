namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * This class is the default envelope-diagram class used by envelope functions
 * to represent the minimization or the maximization diagram of a set of curves
 * in the plane.  It represents the diagram as a doubly-linked list of
 * interleaved vertices and edges. Thus, all operations provided by the envelope
 * diagram take constant time, and the space needed to store the diagram class
 * is linear in the complexity of the envelope.
 *
 * The envelope-diagram class is parameterized by a traits class, which is a
 * model of the `AosXMonotoneTraits_2` concept, in case we handle only envelopes
 * of \f$x\f$-monotone curves, or of the refined `AosTraits_2` concept in case
 * we handle arbitrary planar curves.
 *
 * \cgalModels{EnvelopeDiagram_1}
 *
 * \sa `CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1`
 * \sa `CGAL::Arrangement_on_curve_1::Open_geom_traits_2_adaptor_1`
 * \sa `CGAL::Arrangement_on_curve_1::Open_topology_traits`
 */
template <typename Traits_2, typename Allocator = CGAL_ALLOCATOR(int)>
class Envelope_diagram_1 :
  public Arrangement_on_curve_1::Arrangement_on_curve_1<
    Arrangement_on_curve_1::Open_geom_traits_2_adaptor_1<Traits_2>,
    Arrangement_on_curve_1::Open_topology_traits<Arrangement_on_curve_1::Open_geom_traits_2_adaptor_1<Traits_2>,
                                                 std::list<typename Traits_2::X_monotone_curve_2>,
                                                 std::list<typename Traits_2::X_monotone_curve_2>,
                                                 false, Allocator>> {
public:

  /// \name Types
  /// @{

  //! The container of the curves associated with a diagram feature.
  typedef std::list<typename Traits_2::X_monotone_curve_2> Curve_container;

  //! The 1D geometry traits.
  typedef Arrangement_on_curve_1::Open_geom_traits_2_adaptor_1<Traits_2> Geom_traits_1;

  //! The topology traits of the base 1D arrangement.
  typedef Arrangement_on_curve_1::Open_topology_traits<Geom_traits_1, Curve_container, Curve_container, false,
                                                       Allocator> Topol_traits;

  //! The type of the base 1D arrangement.
  typedef Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits_1, Topol_traits> Arrangement_on_curve_1;

  /// @}

  /// \name Creation
  /// @{

  /*! constructs an empty diagram that shares the 1D geometry traits `traits`.
   * Unlike the constructor from a shared pointer to the 2D geometry traits, it
   * does not allocate a 1D geometry traits.
   */
  Envelope_diagram_1(std::shared_ptr<const Geom_traits_1> traits);

  /// @}

}; // end Envelope_diagram_1

} // end namespace Envelope_2
} // end namespace CGAL
