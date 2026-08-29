// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \ingroup PkgArrangementOnCurve1ConceptsTopologyTraits
 * \cgalConcept
 *
 * The concept `AocTopologyTraits` encapsulates the topological data structures,
 * storage layout, memory allocation, and navigation primitives required by
 * `CGAL::Arrangement_on_curve_1` to maintain subdivisions of a 1D carrier curve.
 *
 * \cgalRefines{DefaultConstructible,CopyConstructible,Assignable}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<Point_1,VertexData,EdgeData,UseVector,Allocator>}
 * \cgalHasModelsEnd
 */
class AocTopologyTraits {
public:
  /// \name Types
  /// @{

  /// The geometric point type stored at vertices.
  typedef unspecified_type Point_1;

  /// Handle type used to refer to a vertex.
  typedef unspecified_type Vertex_descriptor;

  /// Handle type used to refer to an edge.
  typedef unspecified_type Edge_descriptor;

  /// Constant handle type used to refer to a vertex.
  typedef unspecified_type Vertex_const_descriptor;

  /// Constant handle type used to refer to an edge.
  typedef unspecified_type Edge_const_descriptor;

  /// An iterator over vertex descriptors.
  typedef unspecified_type Vertex_descriptor_iterator;

  /// An iterator over edge descriptors.
  typedef unspecified_type Edge_descriptor_iterator;

  /// A range over all vertex descriptors in the arrangement.
  typedef unspecified_type Vertex_descriptor_range;

  /// A range over all edge descriptors in the arrangement.
  typedef unspecified_type Edge_descriptor_range;

  /// Size type.
  typedef unspecified_type Size;

  /// Allocator type.
  typedef unspecified_type Allocator_type;

  /// A readable and writable lvalue property map mapping `Vertex_descriptor` to `const Point_1&`.
  typedef unspecified_type Vertex_point_map;

  /// A readable and writable property map mapping `Vertex_descriptor` to extended vertex user data.
  typedef unspecified_type Vertex_data_map;

  /// A readable and writable property map mapping `Edge_descriptor` to extended edge user data.
  typedef unspecified_type Edge_data_map;

  /// @}

  /// \name Constants
  /// @{

  /// Indicates whether vertices and edges are stored in random-access vector structures.
  static const bool use_vector;

  /// @}

  /// \name Creation
  /// @{

  /// Default constructor. Constructs an empty topology traits.
  AocTopologyTraits();

  /// Constructs an empty topology traits using the specified allocator.
  explicit AocTopologyTraits(const Allocator_type& alloc);

  /// @}

  /// \name Queries and Ranges
  /// @{

  /// returns whether the arrangement contains no vertices.
  bool empty() const;

  /// returns the number of vertices in the arrangement.
  Size number_of_vertices() const;

  /// returns the number of edges in the arrangement.
  Size number_of_edges() const;

  /// returns a range of all vertex descriptors.
  Vertex_descriptor_range vertices() const;

  /// returns a range of all edge descriptors.
  Edge_descriptor_range edges() const;

  /// obtains the property map for vertex geometric points.
  Vertex_point_map vertex_point_map() const;

  /// obtains the property map for vertex user data.
  Vertex_data_map vertex_data_map() const;

  /// obtains the property map for edge user data.
  Edge_data_map edge_data_map() const;

  /// returns the allocator used by the traits.
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

  /// returns whether edge `e` has a left endpoint (is bounded on the left).
  bool has_left_vertex(Edge_const_descriptor e) const;

  /// returns whether edge `e` has a right endpoint (is bounded on the right).
  bool has_right_vertex(Edge_const_descriptor e) const;

  /// directly returns the geometric point associated with vertex `v`.
  const Point_1& vertex_point(Vertex_const_descriptor v) const;

  /// @}

  /// \name Low-Level Storage Primitives
  /// @{

  /// creates a new vertex storing geometric point `p`.
  Vertex_descriptor create_vertex(const Point_1& p);

  /// creates a new edge.
  Edge_descriptor create_edge();

  /// sets `e` as the left incident edge of vertex `v`.
  void set_left_edge(Vertex_descriptor v, Edge_descriptor e);

  /// sets `e` as the right incident edge of vertex `v`.
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e);

  /// sets `v` as the left endpoint vertex of edge `e`.
  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v);

  /// sets `v` as the right endpoint vertex of edge `e`.
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v);

  /// sets `e` as the leftmost unbounded edge.
  void set_unbounded_left_edge(Edge_descriptor e);

  /// sets `e` as the rightmost unbounded edge.
  void set_unbounded_right_edge(Edge_descriptor e);

  /// clears the left vertex endpoint of edge `e`.
  void clear_left_vertex(Edge_descriptor e);

  /// clears the right vertex endpoint of edge `e`.
  void clear_right_vertex(Edge_descriptor e);

  /// erases vertex `v` from the topology structure.
  void erase_vertex(Vertex_descriptor v);

  /// erases edge `e` from the topology structure.
  void erase_edge(Edge_descriptor e);

  /// removes all vertices and resets edges to a single unbounded edge.
  void clear();

  /// swaps the content of the topology traits with `other`.
  void swap(AocTopologyTraits& other) noexcept;

  /// @}
};
