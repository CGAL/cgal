// Copyright (c) 2023 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori

#ifndef CGAL_TRIANGULATION_WITH_EVEN_ODD_CONSTRAINTS_2_H
#define CGAL_TRIANGULATION_WITH_EVEN_ODD_CONSTRAINTS_2_H

#include <CGAL/license/Polygon_repair.h>
#include <CGAL/iterator.h>

namespace CGAL {
namespace Polygon_repair {
namespace internal {

template <class Triangulation_>
class Triangulation_with_even_odd_constraints_2 : public Triangulation_ {
public:
  using Base_triangulation = Triangulation_;
  using Point = typename Triangulation_::Point;
  using Edge = typename Triangulation_::Edge;
  using Vertex_handle = typename Triangulation_::Vertex_handle;
  using Face_handle = typename Triangulation_::Face_handle;
  using List_edges = typename Triangulation_::List_edges;
  using List_faces = typename Triangulation_::List_faces;
  using All_faces_iterator = typename Triangulation_::All_faces_iterator;

  class Interior_tester {
    const Triangulation_with_even_odd_constraints_2 *t;
  public:
    Interior_tester() {}
    Interior_tester(const Triangulation_with_even_odd_constraints_2 *tr) : t(tr) {}

    bool operator()(const All_faces_iterator & fit) const {
      return fit->label() < 1;
    }
  };

  // iterator over interior faces.
  class Interior_faces_iterator : public Filter_iterator<All_faces_iterator, Interior_tester> {
    using Base = Filter_iterator<All_faces_iterator, Interior_tester>;
    using Self = Interior_faces_iterator;
  public:
    Interior_faces_iterator() : Base() {}
    Interior_faces_iterator(const Base &b) : Base(b) {}
    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator Face_handle() const { return Base::base(); }
  };

  // Inserts point p in the triangulation and returns the corresponding vertex.
  Vertex_handle insert(const Point &p, Face_handle f = Face_handle()) {
    return Base_triangulation::insert(p, f);
  }

  // Add constraint from va to vb using the odd-even rule
  void even_odd_insert_constraint(Vertex_handle va, Vertex_handle vb) {

    // Degenerate edge
    if (va == vb) return;

    // [va, vb] is either an existing edge OR
    // there's an existing shorter edge from va in the direction of vb
    Vertex_handle vc; // [va, vc] is the first edge along [va, vb]
    Face_handle incident_face; // incident to [va, vc]
    int opposite_vertex; // opposite to [va, vc]
    if (Base_triangulation::includes_edge(va, vb, vc, incident_face, opposite_vertex)) {
      if (Base_triangulation::is_constrained(Edge(incident_face, opposite_vertex))) {
        Base_triangulation::remove_constrained_edge(incident_face, opposite_vertex);
      } else Base_triangulation::mark_constraint(incident_face, opposite_vertex);
      if (vc != vb) even_odd_insert_constraint(vc, vb); // process edges along [vc, vb]
      return;
    }

    // [va, vb] intersects a constrained edge or an existing vertex
    List_faces intersected_faces;
    List_edges conflict_boundary_ab, conflict_boundary_ba;
    Vertex_handle intersection;
    if (Base_triangulation::find_intersected_faces(va, vb, intersected_faces, conflict_boundary_ab, conflict_boundary_ba, intersection)) {
      if (intersection != va && intersection != vb) {
        even_odd_insert_constraint(va, intersection);
        even_odd_insert_constraint(intersection, vb);
      } else even_odd_insert_constraint(va, vb);
      return;
    }

    // Otherwise
    Base_triangulation::triangulate_hole(intersected_faces, conflict_boundary_ab, conflict_boundary_ba);
    if (intersection != vb) {
      even_odd_insert_constraint(intersection, vb);
    }
  }

  // Add constraint from pa to pb using the odd-even rule
  void even_odd_insert_constraint(Point pa, Point pb) {
    Vertex_handle va = insert(pa);
    Vertex_handle vb = insert(pb, va->face()); // vb is likely close to va
    even_odd_insert_constraint(va, vb);
  }

  // Starts at an arbitrary interior face
  Interior_faces_iterator interior_faces_begin() {
    return CGAL::filter_iterator(Base_triangulation::all_faces_end(),
                                 Interior_tester(this),
                                 Base_triangulation::all_faces_begin());
  }

  // Past-the-end iterator
  Interior_faces_iterator interior_faces_end() {
    return CGAL::filter_iterator(Base_triangulation::all_faces_end(),
                                 Interior_tester(this));
  }
};

} // namespace internal
} // namespace Polygon_repair
} // namespace CGAL

#endif // CGAL_TRIANGULATION_WITH_EVEN_ODD_CONSTRAINTS_2_H
