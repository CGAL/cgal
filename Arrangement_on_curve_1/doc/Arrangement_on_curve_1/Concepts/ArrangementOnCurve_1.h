// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \ingroup PkgArrangementOnCurve1Concepts
 * \cgalConcept
 *
 * \cgalHasModelsBegin
 * \cgalRefines{DefaultConstructible}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>}
 * \cgalHasModelsEnd
 *
 * A model of the concept `ArrangementOnCurve_1` can be used to represent a 1D
 * subdivision of a continuous geometric curve (a 1D "master curve" or line
 * space) into alternating vertices and edges, collectively referred to as cells.
 */
class ArrangementOnCurve_1 {
public:
  /// \name Types
  /// @{

  //! the size type (convertible to `size_t`).
  typedef unspecified_type Size;

  //! The 1D arrangement point representation.
  typedef unspecified_type Point_1;

  //! Descriptor targeting a mutable vertex element.
  typedef unspecified_type Vertex_descriptor;

  //! Descriptor targeting a mutable edge element.
  typedef unspecified_type Edge_descriptor;

  //! Descriptor targeting an immutable vertex element.
  typedef unspecified_type Vertex_const_descriptor;

  //! Descriptor targeting an immutable edge element.
  typedef unspecified_type Edge_const_descriptor;

  /*! Represents a point localization result within the 1D curve structure.
   * A coordinate query point must match a `Vertex_descriptor` if it rests on an arrangement node,
   * or an `Edge_descriptor` if it falls securely between two nodes or out inside an infinity frontier.
   */
  typedef std::variant<Vertex_descriptor, Edge_descriptor> Location_result;
  typedef std::variant<Vertex_const_descriptor, Edge_const_descriptor> Const_location_result;
  /// @}

  /// \name Accessors
  /// @{

  /*! obtains `true` if the arrangement does not contains any vertices.
   */
  bool empty() const;

  /*! obtains the total count of vertices resting along the line track.
   */
  Size number_of_vertices() const;

  /*! obtains the total number of bounded and unbounded edge sections.
   */
  Size number_of_edges() const;

  /*! obtains an iterator range tracking all constant vertex descriptors.
   */
  auto vertices() const;

  /*! obtains an iterator range tracking all constant edge descriptors.
   */
  auto edges() const;

  /*! obtains an lvalue property map matching vertex descriptors to their coordinate positions.
   */
  auto vertex_point_map() const;

  /*! obtains an lvalue property map matching vertex descriptors to their extended user attributes.
   */
  auto vertex_data_map() const;

  /*! obtains an lvalue property map matching edge descriptors to their extended user attributes.
   */
  auto edge_data_map() const;

  /*! obtains the descriptor of the leftmost unbounded edge spanning \f$(-\infty, v_{first})\f$.
   */
  Edge_descriptor unbounded_left_edge();

  /*! obtains the constant descriptor of the leftmost unbounded edge.
   */
  Edge_const_descriptor unbounded_left_edge() const;

  /*! obtains the descriptor of the rightmost unbounded edge spanning \f$(v_{last}, +\infty)\f$.
   */
  Edge_descriptor unbounded_right_edge();

  /*! obtains the constant descriptor of the rightmost unbounded edge.
   */
  Edge_const_descriptor unbounded_right_edge() const;

  /*! obtains the descriptor of the left vertex of an edge.
   * \param e the edge.
   */
  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const { return m_topology_traits.left_vertex(e); }

  /*! obtains the descriptor of the vertex of an edge.
   * \param e the edge.
   */
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return m_topology_traits.right_vertex(e); }

  /*! obtains the descriptor of left edge of a vertex.
   * \param v the vertex.
   * \pre `has_left_vertex(e)` evaluates to `true`.
   */
  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const{ return m_topology_traits.left_edge(v); }

  /*! obtains the descriptor of right edge of a vertex.
   * \param v the vertex.
   * \pre `has_right_vertex(e)` evaluates to `true`.
   */
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const { return m_topology_traits.right_edge(v); }

  /*! determines whether an edge has a left vertex.
   * \param e the edge.
   */
  bool has_left_vertex(Edge_const_descriptor e) const { return m_topology_traits.has_left_vertex(e); }

  /*! determines whether an edge has a right vertex.
   * \param e the edge.
   */
  bool has_right_vertex(Edge_const_descriptor e) const { return m_topology_traits.has_right_vertex(e); }

  /// @}

  /// \name Modifiers
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

  /* inserts the very first point into an empty arrangement.
   * The single unbounded edge is split into (\f$(-\infty\f$, `v`) and (`v`, \f$(+\infty\f$).
   */
  Vertex_descriptor insert_empty(const Point_1& p);

  /*! inserts a new point `p` strictly to the left of an existing vertex `v`.
   */
  Vertex_descriptor insert_before(Vertex_descriptor v, const Point_1& p);

  /*! inserts a new point `p` strictly to the right of an existing vertex `v`.
   */
  Vertex_descriptor insert_after(Vertex_descriptor v, const Point_1& p);

  /*! splits the edge `e` by inserting a new point `p` inside it.
   * `e` becomes the left sub-edge; a new edge becomes the right sub-edge.
   */
  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p);

  /*! removes an active vertex node from the 1D track, cleanly merging its left
   * and right segments into an individual unified edge.
   */
  void remove(Vertex_descriptor v);

  /// @}
};
