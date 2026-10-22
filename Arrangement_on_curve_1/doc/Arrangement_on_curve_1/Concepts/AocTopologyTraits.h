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
 * \cgalRefines{DefaultConstructible,CopyConstructible,Assignable,MoveConstructible,MoveAssignable}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Open_topology_traits<GeometryTraits1, VertexData, EdgeData, UseVector, Allocator>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Closed_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Left_closed_topology_traits<GeometryTraits_1, VertexData, EdgeData,UseVector, Allocator>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Right_closed_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Identified_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>}
 * \cgalHasModelsEnd
 */
class AocTopologyTraits {
public:
  //! \name Types
  //! @{

  //! The allocator type of the internal containers.
  typedef unspecified_type Allocator;

  //! The standard stl type.
  typedef Allocator allocator_type;

  //! Size type.
  typedef std::size_t Size;

  //! The geometric point type stored at vertices.
  typedef unspecified_type Point_1;

  //! The vertex descriptor type (`std::size_t` if `UseVector` is `true`, iterator otherwise).
  typedef unspecified_type Vertex_descriptor;

  //! he edge descriptor type (`std::size_t` if `UseVector` is `true`, iterator otherwise).
  typedef unspecified_type Edge_descriptor;

  //! Constant vertex descriptor type.
  typedef unspecified_type Vertex_const_descriptor;

  //! Constant edge descriptor type.
  typedef unspecified_type Edge_const_descriptor;

  //! An iterator over vertex descriptors.
  typedef unspecified_type Vertex_descriptor_iterator;

  //! An iterator over edge descriptors.
  typedef unspecified_type Edge_descriptor_iterator;

  //! A range over all vertex descriptors in the arrangement.
  typedef unspecified_type Vertex_descriptor_range;

  //! A range over all edge descriptors in the arrangement.
  typedef unspecified_type Edge_descriptor_range;

  //! A readable and writable lvalue property map mapping `Vertex_descriptor` to `const Point_1&`.
  typedef unspecified_type Vertex_point_map;

  //! A readable and writable property map mapping `Vertex_descriptor` to extended vertex user data.
  typedef unspecified_type Vertex_data_map;

  //! A readable and writable property map mapping `Edge_descriptor` to extended edge user data.
  typedef unspecified_type Edge_data_map;

  //! @}

  //! \name Constants
  //! @{

  //! Reflects the value of template parameter `UseVector`.
  static constexpr bool use_vector = UseVector;

  //! @}

  //! \name Queries and Ranges
  //! @{

  //! returns whether the arrangement contains no vertices.
  bool empty() const;

  //! returns the number of vertices in the arrangement.
  Size number_of_vertices() const;

  //! returns the number of edges in the arrangement.
  Size number_of_edges() const;

  //! returns a range of all vertex descriptors.
  Vertex_descriptor_range vertices() const;

  //! returns a range of all edge descriptors.
  Edge_descriptor_range edges() const;

  //! obtains the property map for vertex geometric points.
  Vertex_point_map vertex_point_map() const;

  //! obtains the property map for vertex user data.
  Vertex_data_map vertex_data_map() const;

  //! obtains the property map for edge user data.
  Edge_data_map edge_data_map() const;

  //! returns the allocator used by the traits.
  Allocator get_allocator() const noexcept;

  //! @}

  //! \name Navigation
  //! @{

  //! obtains the leftmost boundary edge.
  Edge_descriptor left_boundary_edge();

  //! obtains the leftmost boundary edge (const version).
  Edge_const_descriptor left_boundary_edge() const;

  //! obtains the rightmost boundary edge.
  Edge_descriptor right_boundary_edge();

  //! obtains the rightmost boundary edge (const version).
  Edge_const_descriptor right_boundary_edge() const;

  //! returns the edge incident to the left of vertex `v`.
  Edge_descriptor left_edge(Vertex_descriptor v);

  //! returns the edge incident to the right of vertex `v`.
  Edge_descriptor right_edge(Vertex_descriptor v);

  //! returns the vertex incident to the left endpoint of edge `e`.
  Vertex_descriptor left_vertex(Edge_descriptor e);

  //! returns the vertex incident to the right endpoint of edge `e`.
  Vertex_descriptor right_vertex(Edge_descriptor e);

  //! returns whether edge `e` has a left endpoint (is closed on the left).
  bool has_left_vertex(Edge_const_descriptor e) const;

  //! returns whether edge `e` has a right endpoint (is closed on the right).
  bool has_right_vertex(Edge_const_descriptor e) const;

  //! directly returns the geometric point associated with vertex `v`.
  const Point_1& vertex_point(Vertex_const_descriptor v) const;

  //! @}

  //! \name Low-Level Storage Primitives
  //! @{

  //! sets the leftmost boundary edge descriptor.
  void set_left_boundary_edge(Edge_descriptor e);

  //! sets the rightmost boundary edge descriptor.
  void set_right_boundary_edge(Edge_descriptor e);

  //! creates a new vertex associated with point `p`.
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt);

  //! creates a new uninitialized edge.
  Edge_descriptor create_edge();

  //! sets `e` as the left incident edge of vertex `v`.
  void set_left_edge(Vertex_descriptor v, Edge_descriptor e);

  //! sets `e` as the right incident edge of vertex `v`.
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e);

  //! sets `v` as the left endpoint vertex of edge `e`.
  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v);

  //! sets `v` as the right endpoint vertex of edge `e`.
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v);

  //! unbinds the left vertex of edge `e`.
  void clear_left_vertex(Edge_descriptor e);

  //! unbinds the right vertex of edge `e`.
  void clear_right_vertex(Edge_descriptor e);

  //! erases vertex `v`.
  void erase_vertex(Vertex_descriptor v);

  /*! erases edge `e`.
   * \pre `e` is neither the left nor the right boundary edge.
   */
  void erase_edge(Edge_descriptor e);

  //! removes all vertices and edges and resets the topology to its initial state.
  void clear();

  //! swaps the content of the topology traits with `other`.
  void swap(AocTopologyTraits& other) noexcept;

  //! @}
};
