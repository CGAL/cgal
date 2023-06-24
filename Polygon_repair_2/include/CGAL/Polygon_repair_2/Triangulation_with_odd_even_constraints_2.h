// Copyright (c) 2023 GeometryFactory. All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori

#ifndef CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H
#define CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H

namespace CGAL {

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Triangulation_with_odd_even_constraints_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 *
 * \tparam Triangulation_ must have support for constraints
 */
template <class Triangulation_>
class Triangulation_with_odd_even_constraints_2 : public Triangulation_ {
public:
  /// \name Definition

  /// @{
  /// the triangulation class.
  typedef Triangulation_ Base_triangulation;

  /// Point type
  typedef typename Triangulation_::Point Point;

  /// Edge type
  typedef typename Triangulation_::Edge Edge;

  /// handle to a vertex
  typedef typename Triangulation_::Vertex_handle Vertex_handle;

  /// handle to a face
  typedef typename Triangulation_::Face_handle Face_handle;

  /// list of edges
  typedef typename Triangulation_::List_edges List_edges;
  
  /// list of faces
  typedef typename Triangulation_::List_faces List_faces;

  /// @}
  
  // iterator over interior faces.
  class Interior_faces_iterator : public Base_triangulation::All_faces_iterator {
    Interior_faces_iterator operator++();
    Interior_faces_iterator operator--();
  };

  // Inserts point p in the triangulation and returns the corresponding vertex.
  Vertex_handle insert(const Point &p, Face_handle f = Face_handle()) {
    return Base_triangulation::insert(p, f);
  }

  // Add constraint from va to vb using the odd-even rule
  void odd_even_insert_constraint(Vertex_handle va, Vertex_handle vb) {

    // Degenerate edge
    if (va == vb) return;

    // [va, vb] is an existing edge OR it is composed of shorter edges
    Vertex_handle vc; // [va, vc] is the first edge along [va, vb]
    Face_handle incident_face; // incident to [va, vc]
    int opposite_vertex; // opposite to [va, vc]
    if (Base_triangulation::includes_edge(va, vb, vc, incident_face, opposite_vertex)) {
      if (Base_triangulation::is_constrained(Edge(incident_face, opposite_vertex))) {
        Base_triangulation::remove_constrained_edge(incident_face, opposite_vertex);
      } else Base_triangulation::mark_constraint(incident_face, opposite_vertex);
      if (vc != vb) odd_even_insert_constraint(vc, vb); // process edges along [vc, vb]
      return;
    }

    // [va, vb] intersects a constrained edge or an existing vertex
    List_faces intersected_faces;
    List_edges conflict_boundary_ab, conflict_boundary_ba;
    Vertex_handle intersection;
    if (Base_triangulation::find_intersected_faces(va, vb, intersected_faces, conflict_boundary_ab, conflict_boundary_ba, intersection)) {
      if (intersection != va && intersection != vb) {
        odd_even_insert_constraint(va, intersection);
        odd_even_insert_constraint(intersection, vb);
      } else odd_even_insert_constraint(va, vb);
      return;
    }

    // Otherwise
    Base_triangulation::triangulate_hole(intersected_faces, conflict_boundary_ab, conflict_boundary_ba);
    if (intersection != vb) {
      odd_even_insert_constraint(intersection, vb);
    }
  }

  // Add constraint from pa to pb using the odd-even rule
  void odd_even_insert_constraint(Point pa, Point pb) {
    Vertex_handle va = insert(pa);
    Vertex_handle vb = insert(pb);
    odd_even_insert_constraint(va, vb);
  }

  // Starts at an arbitrary interior face
  Interior_faces_iterator interior_faces_begin() {

  }

  // Past-the-end iterator
  Interior_faces_iterator interior_faces_end() {
    
  }
};

} // namespace CGAL

#endif // CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H