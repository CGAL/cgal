// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Ref
 *
 * The class template `Unbounded_topology_traits` provides a model of the
 * `AocTopologyTraits` concept for 1D arrangements embedded on an unbounded
 * carrier curve.
 *
 * It maintains the arrangement as an alternating sequence of vertices and edges
 * flanked by two semi-infinite unbounded edges (`unbounded_left_edge` and
 * `unbounded_right_edge`). The class template supports both list-based and
 * vector-based underlying storage structures, and permits optional user data
 * extension on vertices and edges via property maps.
 *
 * \cgalModels{AocTopologyTraits}
 *
 * \tparam Point_1 the geometric point type stored at vertices.
 * \tparam VertexData optional user-defined data attached to each vertex (defaults to `void`).
 * \tparam EdgeData optional user-defined data attached to each edge (defaults to `void`).
 * \tparam UseVector when set to `true`, stores vertices and edges in `std::vector`
 *   to enable \f$O(\log n)\f$ binary search point location. Descriptors become `std::size_t`
 *   indices. When `false` (default), uses `std::list` where descriptors are stable
 *   bidirectional iterators.
 * \tparam Allocator allocator used for internal node and container allocations (defaults to `std::allocator<char>`).
 *
 * \sa `AocTopologyTraits`
 * \sa `CGAL::Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>`
 */
template <typename Point_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
class Unbounded_topology_traits {
public:
  /// \name Types
  /// @{

  /// Allocator type.
  typedef Allocator Allocator_type;

  /// Size type.
  typedef std::size_t Size;

  /// The geometric point type.
  typedef Point_1 Point_1;

  /// The vertex descriptor type (`std::size_t` if `UseVector` is `true`, iterator otherwise).
  typedef unspecified_type Vertex_descriptor;

  /// The edge descriptor type (`std::size_t` if `UseVector` is `true`, iterator otherwise).
  typedef unspecified_type Edge_descriptor;

  /// Constant vertex descriptor type.
  typedef unspecified_type Vertex_const_descriptor;

  /// Constant edge descriptor type.
  typedef unspecified_type Edge_const_descriptor;

  /// Iterator over vertex descriptors.
  typedef unspecified_type Vertex_descriptor_iterator;

  /// Iterator over edge descriptors.
  typedef unspecified_type Edge_descriptor_iterator;

  /// Range of vertex descriptors.
  typedef boost::iterator_range<Vertex_descriptor_iterator> Vertex_descriptor_range;

  /// Range of edge descriptors.
  typedef boost::iterator_range<Edge_descriptor_iterator> Edge_descriptor_range;

  /// Lvalue property map mapping `Vertex_descriptor` to `const Point_1&`.
  typedef unspecified_type Vertex_point_map;

  /// Lvalue property map mapping `Vertex_descriptor` to `VertexData&`.
  typedef unspecified_type Vertex_data_map;

  /// Lvalue property map mapping `Edge_descriptor` to `EdgeData&`.
  typedef unspecified_type Edge_data_map;

  /// @}

  /// \name Constants
  /// @{

  /// Reflects the value of template parameter `UseVector`.
  static constexpr bool use_vector = UseVector;

  /// @}

  /// \name Creation
  /// @{

  /// constructs an empty topology traits initialized with a single unbounded edge.
  Unbounded_topology_traits();

  /// constructs an empty topology traits with an explicit allocator instance.
  explicit Unbounded_topology_traits(const Allocator_type& alloc);

  /// copy constructor.
  Unbounded_topology_traits(const Unbounded_topology_traits& other);

  /// move constructor.
  Unbounded_topology_traits(Unbounded_topology_traits&& other) noexcept;

  /// @}

  /// \name Queries and Ranges
  /// @{

  /// returns `true` if the arrangement contains no vertices, and `false` otherwise.
  bool empty() const;

  /// returns the number of vertices in the arrangement.
  Size number_of_vertices() const;

  /// returns the number of edges in the arrangement.
  Size number_of_edges() const;

  /// returns an iterator range yielding all vertex descriptors.
  Vertex_descriptor_range vertices() const;

  /// returns an iterator range yielding all edge descriptors.
  Edge_descriptor_range edges() const;

  /// obtains the property map for accessing vertex geometric coordinates.
  Vertex_point_map vertex_point_map() const;

  /// obtains the property map for accessing vertex user data.
  Vertex_data_map vertex_data_map() const;

  /// obtains the property map for accessing edge user data.
  Edge_data_map edge_data_map() const;

  /// returns the allocator instance.
  Allocator_type get_allocator() const noexcept;

  /// @}

  /// \name Navigation
  /// @{

  /// obtains the leftmost unbounded edge.
  Edge_descriptor unbounded_left_edge();

  /// obtains the leftmost unbounded edge (const version).
  Edge_const_descriptor unbounded_left_edge() const;

  /// obtains the rightmost unbounded edge.
  Edge_descriptor unbounded_right_edge();

  /// obtains the rightmost unbounded edge (const version).
  Edge_const_descriptor unbounded_right_edge() const;

  /// returns the edge incident to the left of vertex `v`.
  Edge_descriptor left_edge(Vertex_descriptor v);

  /// returns the edge incident to the right of vertex `v`.
  Edge_descriptor right_edge(Vertex_descriptor v);

  /// returns the vertex incident to the left endpoint of edge `e`.
  Vertex_descriptor left_vertex(Edge_descriptor e);

  /// returns the vertex incident to the right endpoint of edge `e`.
  Vertex_descriptor right_vertex(Edge_descriptor e);

  /// returns whether edge `e` has a left endpoint.
  bool has_left_vertex(Edge_const_descriptor e) const;

  /// returns whether edge `e` has a right endpoint.
  bool has_right_vertex(Edge_const_descriptor e) const;

  /// returns the geometric point associated with vertex `v` without property map overhead.
  const Point_1& vertex_point(Vertex_const_descriptor v) const;

  /// @}

  /// \name Modifiers
  /// @{

  /// creates a new vertex associated with point `p` and appends it to the container.
  Vertex_descriptor create_vertex(const Point_1& p);

  /// creates a new uninitialized edge.
  Edge_descriptor create_edge();

  /// sets `e` as the left incident edge of vertex `v`.
  void set_left_edge(Vertex_descriptor v, Edge_descriptor e);

  /// sets `e` as the right incident edge of vertex `v`.
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e);

  /// sets `v` as the left endpoint of edge `e`.
  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v);

  /// sets `v` as the right endpoint of edge `e`.
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v);

  /// sets the leftmost unbounded edge descriptor.
  void set_unbounded_left_edge(Edge_descriptor e);

  /// sets the rightmost unbounded edge descriptor.
  void set_unbounded_right_edge(Edge_descriptor e);

  /// unbinds the left vertex of edge `e`.
  void clear_left_vertex(Edge_descriptor e);

  /// unbinds the right vertex of edge `e`.
  void clear_right_vertex(Edge_descriptor e);

  /// erases vertex `v`.
  void erase_vertex(Vertex_descriptor v);

  /// erases edge `e`.
  void erase_edge(Edge_descriptor e);

  /// clears all vertices and resets edges to a single unbounded edge.
  void clear();

  /// swaps with `other`.
  void swap(Unbounded_topology_traits& other) noexcept;

  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
