// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_OVERLAY_VISITOR_H
#define CGAL_OVERLAY_VISITOR_H

namespace CGAL {

/*! \ingroup PkgArrangementOnCurve1Concepts
 * \cgalConcept
 *
 * A model of the `OverlayVisitor` concept contains callback functions invoked
 * while the overlay operation of two input arrangements progresses.  Models of
 * the concept are used by the free function `Arrangement_on_curve_1::overlay()`
 * to maintain the auxiliary data stored with the cells (i.e., vertices and
 * edges) of the resulting overlaid arrangement, based on the contents of the
 * input cells.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Default_overlay_visitor<ArrangementA, ArrangementB, ArrangementR>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Copy_overlay_visitor<ArrangementA, ArrangementB, ArrangementR>}
 * \cgalHasModelsEnd
 *
 * \sa `overlay`
 */
class OverlayVisitor {
public:
  //! a constant descriptor of a vertex in the first input arrangement.
  typedef Vertex_const_descriptor_a;

  //! a constant descriptor of a vertex in the second input arrangement.
  typedef Vertex_const_descriptor_b;

  //! a descriptor of a vertex in the output arrangement.
  typedef Vertex_descriptor_r;

  //! a constant descriptor of an edge in the first input arrangement.
  typedef Edge_const_descriptor_a;

  //! a constant descriptor of an edge in the second input arrangement.
  typedef Edge_const_descriptor_b;

  //! a descriptor of an edge in the output arrangement.
  typedef Edge_descriptor_r;

  /*! updates the vertex identified by the descriptor `v_res`, which is induced
   * by the vertices of the two input arrangements identified by the descriptors
   * `v_a` and `v_b`, repsectively.
   */
  void create_vertex(Vertex_const_descriptor_a v_a, Vertex_const_descriptor_b v_b, Vertex_descriptor_r v_res);

  /*! updates the vertex identified by the descriptor `v_res`, which is induced
   * by tyhe vertex of the first input arrangements and the edge of the second
   * input arrangement, identified by the descriptors `v_a` and `e_b`,
   * repsectively.
   */
  void create_vertex(Vertex_const_descriptor_a v_a, Edge_const_descriptor_b e_b, Vertex_descriptor_r v_res);

  /*! updates the vertex identified by the descriptor `v_res`, which is induced
   * by the edge of the first input arrangements and the vertex of the second
   * input arrangement, identified by the descriptors `e_a` and `v_b`,
   * repsectively.
   */
  void create_vertex(Edge_const_descriptor_a e_a, Vertex_const_descriptor_b v_b, Vertex_descriptor_r v_res);

  /*! updates the vertex identified by the descriptor `v_res`, which is induced
   * by the edges of the two input arrangements identified by the descriptors
   * `e_a` and `e_b`, repsectively.
   */
  void create_edge(Edge_const_descriptor_a e_a, Edge_const_descriptor_b e_b, Edge_descriptor_r e_res);
}

} // namespace CGAL

#endif
