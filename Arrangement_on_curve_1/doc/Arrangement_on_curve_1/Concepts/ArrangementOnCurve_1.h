// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \ingroup PkgArrangementOnCurve1Concepts
 * \cgalConcept
 *
 * \cgalRefines{DefaultConstructible}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>}
 * \cgalHasModelsEnd
 *
 * A model of the concept `ArrangementOnCurve_1` can be used to represent a 1D
 * subdivision of a continuous geometric curve (a 1D "master curve" or line
 * space) into alternating vertices and edges, collectively referred to as cells.
 *
 * The geometry is governed by the `GeometryTraits_1` template parameter, which
 * must model the concept `AocTraits_1`.  In particular, any model of the concept
 * `GeometryTraits_1` must define the `Left_side_category` and
 * `Right_side_category` categories to describe whether the parameter space is
 * open, closed, or identified at each end.  When both categories are
 * `CGAL::Arrangement_on_curve_1::Identified_side_tag`, the curve is cyclic.
 */
class ArrangementOnCurve_1 {
public:
  //! \name Types
  //! @{

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
   * A coordinate query point matches a `Vertex_descriptor` if it rests on an
   * arrangement vertex, or an `Edge_descriptor` if it falls inside an edge.
   */
  typedef std::variant<Vertex_descriptor, Edge_descriptor> Location_result;
  typedef std::variant<Vertex_const_descriptor, Edge_const_descriptor> Const_location_result;

  //! @}

  //! \name Boundary Categories
  //! @{

  //! Tag type describing the left boundary condition.
  //! Reflects `GeometryTraits_1::Left_side_category`.
  typedef unspecified_type Left_side_category;

  //! Tag type describing the right boundary condition.
  //! Reflects `GeometryTraits_1::Right_side_category`.
  typedef unspecified_type Right_side_category;

  //! @}

  //! \name Accessors
  //! @{

  /*! is `true` if the curve is cyclic (both
   * `Left_side_category` and `Right_side_category` are `CGAL::Arrangement_on_curve_1::Identified_side_tag`).
   */
  static constexpr bool is_identified;

  /*! returns `true` if the arrangement does not contain any vertices.
   */
  bool empty() const;

  /*! obtains the total count of vertices resting along the carrier curve.
   */
  Size number_of_vertices() const;

  /*! obtains the total number of edges.
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

  /*! obtains the descriptor of the leftmost edge.  If the parameter space is
   * identified, this is the edge whose left vertex is the identification
   * vertex, or, if the identification vertex does not exist, the edge that
   * contains the identification point.
   */
  Edge_descriptor left_boundary_edge();

  /*! obtains the constant descriptor of the leftmost edge.
   */
  Edge_const_descriptor left_boundary_edge() const;

  /*! obtains the descriptor of the rightmost edge.  If the parameter space is
   * identified, this is the edge whose right vertex is the identification
   * vertex, or, if the identification vertex does not exist, the edge that
   * contains the identification point, that is,
   * \link ArrangementOnCurve_1::left_boundary_edge() `left_boundary_edge()`\endlink.
   */
  Edge_descriptor right_boundary_edge();

  /*! obtains the constant descriptor of the rightmost edge.
   */
  Edge_const_descriptor right_boundary_edge() const;

  /*! obtains the descriptor of the left vertex of an edge.
   * \param e the edge.
   * \pre `has_left_vertex(e)` is `true`.
   */
  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const;

  /*! obtains the descriptor of the right vertex of an edge.
   * \param e the edge.
   * \pre `has_right_vertex(e)` is `true`.
   */
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const;

  /*! obtains the descriptor of the left edge of a vertex.
   * \param v the vertex.
   */
  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const;

  /*! obtains the descriptor of the right edge of a vertex.
   * \param v the vertex.
   */
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const;

  /*! determines whether an edge has a left vertex.
   * \param e the edge.
   */
  bool has_left_vertex(Edge_const_descriptor e) const;

  /*! determines whether an edge has a right vertex.
   * \param e the edge.
   */
  bool has_right_vertex(Edge_const_descriptor e) const;

  //! @}

  //! \name Modifiers
  //! @{

  /*! creates a new vertex, enforcing the sorted ordering invariant when
   * `BinarySearch` is `true`.
   */
  Vertex_descriptor create_vertex(const Point_1& p);

  /*! creates a new edge.
   */
  Edge_descriptor create_edge();

  /*! destroys a given vertex.
   */
  void destroy_vertex(Vertex_descriptor v);

  /*! destroys a given edge.
   */
  void destroy_edge(Edge_descriptor e);

  /*! inserts the very first point into an empty arrangement.  If the
   * parameter space is not identified, the single edge that spans the
   * parameter space is split into the edge to the left of the new vertex `v`
   * and the edge to the right of `v`.  Otherwise, the single edge becomes a
   * self-loop whose two ends are `v`.
   *
   * \param p the point to insert.
   */
  Vertex_descriptor insert_empty(const Point_1& p);

  /*! inserts a new point `p` strictly to the left of an existing vertex `v`.
   * \pre `p` lies in the interior of the parameter space; that is, it lies on
   * no closed boundary, and it is not the identification point.
   */
  Vertex_descriptor insert_before(Vertex_descriptor v, const Point_1& p);

  /*! inserts a new point `p` strictly to the right of an existing vertex `v`.
   * \pre `p` lies in the interior of the parameter space; that is, it lies on
   * no closed boundary, and it is not the identification point.
   */
  Vertex_descriptor insert_after(Vertex_descriptor v, const Point_1& p);

  /*! splits the edge `e` by inserting a new point `p` inside it.
   * `e` becomes the left sub-edge; a new edge becomes the right sub-edge.
   * \pre `p` lies in the interior of the parameter space; that is, it lies on
   * no closed boundary, and it is not the identification point.
   */
  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p);

  /*! removes vertex `v` from the arrangement.  If the two edges incident to
   * `v` are distinct, they are merged into a single edge.  Otherwise, `v` is
   * a boundary vertex of a closed side, or the only vertex of an arrangement
   * on an identified parameter space, and no edge is removed; in the latter
   * case the arrangement becomes empty.
   */
  void remove(Vertex_descriptor v);

  /*! clears the arrangement.
   */
  void clear();

  /*! reserves storage for `num_vertices` vertices and the edges between them,
   * so that, in vector mode, no reallocation occurs until the arrangement has
   * more vertices; has no effect in list mode.
   */
  void reserve(Size num_vertices);

  //! @}
};
