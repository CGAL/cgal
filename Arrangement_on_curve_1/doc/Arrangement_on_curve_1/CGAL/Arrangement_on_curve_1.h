namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Ref
 *
 * \cgalModels{ArrangementOnCurve_1}
 *
 * \tparam GeometryTraits_1 must be a model of the `AocTraits_1` concept, providing the geometric
 *         types and a `Compare_x_1` functor to sort points linearly along the continuous curve trajectory.
 * \tparam TopologyTraits must be a model of the 1D arrangement topology traits concept,
 *         managing container allocations, structural adjacency records, and property mappings.
 * \tparam BinarySearch when true (requires `TopologyTraits::UseVector = true`), the `locate()` function
 *         uses a binary search over the sorted vertex vector, reducing locate from \f$O(n)\f$ to \f$O(\log n)\f$.
 *         However, insertion and removal of vertices is limited to rightmost vertices only.
 *
 * An object, the type of which is an instance of the class template
 * `Arrangement_on_curve_1<GeometryTraits_1,TopologyTraits,BinarySearch>`,
 * represents the subdivision as a linked list of alternating vertices and
 * edges.  Thus, the operations that obtain incident cells provided by the class
 * template take constant time, and the space needed to store the arrangement is
 * linear in the complexity of the arrangement.
 *
 * An arrangement object maintains a smart pointer to a shared immutable
 * geometry traits object, optimizing lifecycle tracking and enabling efficient
 * data synchronization across topological overlays.
 *
 * \sa `AocTraits_1`
 */
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch = false>
class Arrangement_on_curve_1 {
public:
  /// \name Types
  /// @{

  //! The geometry traits type.
  typedef GeometryTraits_1 Geometry_traits_1;

  //! The topology traits type.

  typedef TopologyTraits Topology_traits;

  //! A smart pointer to the shared geometry traits context.
  typedef std::shared_ptr<const Geometry_traits_1> Shared_geometry_traits;

  /// @}

  /// \name Creation
  /// @{

  /*! Constructor from an existing geometry traits smart pointer.
   * \param shared_geometry_traits A shared pointer to an existing, valid geometry traits instance.
   */
  explicit Arrangement_on_curve_1(const Shared_geometry_traits shared_geometry_traits);

  /*! Fully custom constructor passing a traits pointer and an explicit topology configuration block.
   */
  Arrangement_on_curve_1(const Shared_geometry_traits shared_geometry_traits, Topology_traits topo_tr);

  /// @}

  /// \name Accessors
  /// @{

  /*! obtains a reference to the active geometry traits object.
   */
  const Geometry_traits_1& geometry_traits_1() const;

  /*! obtains the shared smart pointer tracking the geometry traits context.
   */
  Shared_geometry_traits shared_geometry_traits_1() const;

  /*! obtains an immutable reference to the internal topology traits object.
   */
  const Topology_traits& topology_traits() const;

  /*! obtains a mutable reference to the internal topology traits object.
   */
  Topology_traits& topology_traits();

  /// @}

  /// \name Modification Modifiers
  /// @{

  /*! creates a new vertex, enforcing the rightmost ordering invariant when BinarySearch is active.
   */
  Vertex_descriptor create_vertex(const Point_1& p);

  /*! creates an new edge.
   */
  Edge_descriptor create_edge();

  /*! destroys a given vertex.
   */
  void destroy_vertex(Vertex_descriptor v);

  /*! destroys a given edge.
   */
  void destroy_edge(Edge_descriptor e);

  /*! inserts the first structural node into an empty subdivision context;
   * splits the singular unbounded sequence edge into an absolute left-infinity
   * edge and an absolute right-infinity edge.
   */
  Vertex_descriptor insert_empty(const Point_1& p);

  /*! inserts a point strictly to the left of the leftmost vertex.
   */
  Vertex_descriptor insert_before(Vertex_descriptor v, const Point_1& p);

  /*! inserts a point strictly to the right of the rightmost vertex.
   */
  Vertex_descriptor insert_after(Vertex_descriptor v, const Point_1& p);

  /*! splits an existing edge interval at the designated coordinate parameter position.
   */
  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p);

  /*! removes an active vertex node from the 1D track, cleanly merging its left
   * and right segments into an individual unified edge.
   */
  void remove(Vertex_descriptor v);

  /// @}
};

} // end namespace Arrangement_on_curve_1
} // end namespace CGAL
