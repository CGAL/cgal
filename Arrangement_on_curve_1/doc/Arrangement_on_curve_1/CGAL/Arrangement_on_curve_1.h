// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/**
 * \ingroup PkgArrangementOnCurve1Ref
 *
 * \brief An instance of the class `Arrangement_on_curve_1` models a 1D subdivision of a
 * continuous geometric curve (a 1D "master curve" or line space) into vertices and edges.
 *
 * \tparam GeometryTraits_1 must be a model of the `AocTraits_1` concept, providing the geometric
 * types and a `Compare_x_1` functor to sort points linearly along the continuous curve trajectory.
 * \tparam TopologyTraits must be a model of the 1D arrangement topology traits concept,
 * managing container allocations, structural adjacency records, and property mappings.
 *
 * An arrangement object maintains a smart pointer to a shared immutable geometry traits instance,
 * optimizing lifecycle tracking and enabling efficient data synchronization across topological overlays.
 *
 * \sa `AocTraits_1`
 */
template <typename GeometryTraits_1, typename TopologyTraits>
class Arrangement_on_curve_1 {
public:
  /// \name Types
  /// @{
  typedef GeometryTraits_1 Geometry_traits_1; ///< The geometry traits type.
  typedef TopologyTraits Topology_traits;     ///< The topology traits type.
  typedef std::shared_ptr<const Geometry_traits_1> Shared_geometry_traits; ///< A smart pointer to the shared geometry traits context.
  typedef typename Geometry_traits_1::Point_1 Point_1; ///< The 1D arrangement point representation.

  typedef typename Topology_traits::Vertex_descriptor Vertex_descriptor; ///< Descriptor handle targeting a mutable vertex element.
  typedef typename Topology_traits::Edge_descriptor Edge_descriptor;     ///< Descriptor handle targeting a mutable edge element.
  typedef typename Topology_traits::Vertex_const_descriptor Vertex_const_descriptor; ///< Descriptor handle targeting an immutable vertex element.
  typedef typename Topology_traits::Edge_const_descriptor Edge_const_descriptor;     ///< Descriptor handle targeting an immutable edge element.

  /**
   * \brief Represents a point localization result within the 1D curve structure.
   * A coordinate query point must match a `Vertex_descriptor` if it rests on an arrangement node,
   * or an `Edge_descriptor` if it falls securely between two nodes or out inside an infinity frontier.
   */
  typedef std::variant<Vertex_descriptor, Edge_descriptor> Location_result;
  typedef std::variant<Vertex_const_descriptor, Edge_const_descriptor> Const_location_result;
  /// @}

  /// \name Creation
  /// @{

  /**
   * \brief Default Constructor. Allocates a new default-constructed geometry traits instance on the heap.
   */
  Arrangement_on_curve_1();

  /**
   * \brief Constructor from an existing geometry traits smart pointer.
   * \param shared_geometry_traits A shared pointer to an existing, valid geometry traits instance.
   */
  explicit Arrangement_on_curve_1(const Shared_geometry_traits shared_geometry_traits);

  /**
   * \brief Fully custom constructor passing a traits pointer and an explicit topology configuration block.
   */
  Arrangement_on_curve_1(const Shared_geometry_traits shared_geometry_traits, Topology_traits topo_tr);
  /// @}

  /// \name Accessors
  /// @{
  const Geometry_traits_1& geometry_traits_1() const;       ///< Returns a reference to the active geometry traits instance.
  Shared_geometry_traits shared_geometry_traits_1() const;  ///< Returns the shared smart pointer tracking the geometry traits context.
  const Topology_traits& topology_traits() const;           ///< Returns an immutable reference to the internal topology store.
  Topology_traits& topology_traits();                       ///< Returns a mutable reference to the internal topology store.

  bool is_empty() const;                   ///< Returns `true` if the arrangement contains no vertices.
  size_t number_of_vertices() const;       ///< Returns the total count of vertices resting along the line track.
  size_t number_of_edges() const;          ///< Returns the total number of bounded and unbounded edge sections.

  auto vertices() const;                   ///< Returns an iterator range tracking all constant vertex descriptors.
  auto edges() const;                      ///< Returns an iterator range tracking all constant edge descriptors.

  auto vertex_point_map() const;           ///< Returns an lvalue property map matching vertex descriptors to their coordinate positions.
  auto vertex_data_map() const;            ///< Returns an lvalue property map matching vertex descriptors to their extended user attributes.
  auto edge_data_map() const;              ///< Returns an lvalue property map matching edge descriptors to their extended user attributes.

  Edge_descriptor unbounded_left_edge();               ///< Returns a handle to the leftmost unbounded edge spanning \f$(-\infty, v_{first})\f$.
  Edge_const_descriptor unbounded_left_edge() const;   ///< Returns a constant handle to the leftmost unbounded edge.
  Edge_descriptor unbounded_right_edge();              ///< Returns a handle to the rightmost unbounded edge spanning \f$(v_{last}, +\infty)\f$.
  Edge_const_descriptor unbounded_right_edge() const;  ///< Returns a constant handle to the rightmost unbounded edge.
  /// @}

  /// \name Modification Modifiers
  /// @{

  /**
   * \brief Safety reset function allowing empty arrangements to bind to a separate existing traits memory frame.
   * \pre `is_empty() == true`
   */
  void reset_shared_geometry_traits(Shared_geometry_traits new_shared_traits_traits);

  /**
   * \brief Inserts the first structural node into an empty subdivision context.
   * Splits the singular unbounded sequence edge into an absolute left-infinity edge and an absolute right-infinity edge.
   */
  Vertex_descriptor insert_empty(const Point_1& p);

  /**
   * \brief Inserts a point strictly to the left of the leftmost vertex.
   */
  Vertex_descriptor insert_before(Vertex_descriptor v, const Point_1& p);

  /**
   * \brief Inserts a point strictly to the right of the rightmost vertex.
   */
  Vertex_descriptor insert_after(Vertex_descriptor v, const Point_1& p);

  /**
   * \brief Splits an existing edge interval at the designated coordinate parameter position.
   */
  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p);

  /**
   * \brief Removes an active vertex node from the 1D track, cleanly merging its left and right segments into an individual unified edge.
   */
  void remove(Vertex_descriptor v);
  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
