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

#ifndef CGAL_POLYGON_REPAIR_2_H
#define CGAL_POLYGON_REPAIR_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Polygon_repair_2/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_face_base_with_repair_info_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_with_odd_even_constraints_2.h>

namespace CGAL {

namespace Polygon_repair_2 {

template <class Kernel, class PolygonContainer>
class Polygon_repair_2;

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon without holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Polygon_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Polygon_with_holes_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a multipolygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& mp) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(mp);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Polygon_repair_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 */
template <class Kernel = CGAL::Exact_predicates_inexact_constructions_kernel,
          class PolygonContainer = std::vector<typename Kernel::Point_2>>
class Polygon_repair_2 {
public:
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vertex_base;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
  typedef CGAL::Triangulation_face_base_with_repair_info_2<Kernel, Face_base> Face_base_with_repair_info;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info> Triangulation_data_structure;
  typedef CGAL::Exact_predicates_tag Tag; // assumed for now
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_Delaunay_triangulation;
  typedef Triangulation_with_odd_even_constraints_2<Constrained_Delaunay_triangulation> Triangulation;

  /// \name Creation
  Polygon_repair_2() {

  }

  /// \name Modifiers
  /// @{

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Polygon_2<Kernel, PolygonContainer>& polygon) {
    for (auto const& edge: polygon.edges()) {
      t.odd_even_insert_constraint(edge.source(), edge.target());
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Polygon_with_holes_2<Kernel, PolygonContainer>& polygon) {
    add_to_triangulation(polygon.outer_boundary());
    for (auto const& hole: polygon.holes()) {
       add_to_triangulation(hole);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& multipolygon) {
    for (auto const& polygon: multipolygon.polygons()) {
      add_to_triangulation(polygon);
    }
  }

  // Label triangles in triangulation as inside or outside the polygon
  void label_triangulation() {
    for (auto const face: t.all_face_handles()) {
      face->processed() = false;
      face->interior() = false;
    } std::list<typename Triangulation::Face_handle> to_check;
    t.infinite_face()->processed() = true;
    to_check.push_back(t.infinite_face());
    while (!to_check.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!to_check.front()->neighbor(neighbour)->processed()) {
          to_check.front()->neighbor(neighbour)->processed() = true;
          if (t.is_constrained(typename Triangulation::Edge(to_check.front(), neighbour))) {
            to_check.front()->neighbor(neighbour)->interior() = !to_check.front()->interior();
            to_check.push_back(to_check.front()->neighbor(neighbour));
          } else {
            to_check.front()->neighbor(neighbour)->interior() = to_check.front()->interior();
            to_check.push_back(to_check.front()->neighbor(neighbour));
          }
        }
      } to_check.pop_front();
    }
  }

  void get_boundary(typename Triangulation::Face_handle face, int edge, std::list<typename Triangulation::Vertex_handle> &out_vertices) {
    // Check clockwise edge
    if (face->neighbor(face->cw(edge))->interior() && !face->neighbor(face->cw(edge))->reconstructed()) {
      face->neighbor(face->cw(edge))->reconstructed() = true;
      std::list<typename Triangulation::Vertex_handle> v1;
      get_boundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face), v1);
      out_vertices.splice(out_vertices.end(), v1);
    }
    
    // Add central vertex
    out_vertices.push_back(face->vertex(edge));
    
    // Check counterclockwise edge
    if (face->neighbor(face->ccw(edge))->interior() && !face->neighbor(face->ccw(edge))->reconstructed()) {
      face->neighbor(face->ccw(edge))->reconstructed() = true;
      std::list<typename Triangulation::Vertex_handle> v2;
      get_boundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face), v2);
      out_vertices.splice(out_vertices.end(), v2);
    }
  }

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
    for (auto const face: t.all_face_handles()) {
      face->reconstructed() = false;
    }

    for (auto const interior_face: t.finite_face_handles()) {
      if (!interior_face->interior() || interior_face->reconstructed()) continue;

      // Get boundary
      std::list<typename Triangulation::Vertex_handle> vertices;
      if (interior_face->neighbor(2)->interior() && !interior_face->neighbor(2)->reconstructed()) {
        interior_face->neighbor(2)->reconstructed() = true;
        std::list<typename Triangulation::Vertex_handle> l2;
        get_boundary(interior_face->neighbor(2), interior_face->neighbor(2)->index(interior_face), l2);
        vertices.splice(vertices.end(), l2);
      } vertices.push_back(interior_face->vertex(0));
      if (interior_face->neighbor(1)->interior() && !interior_face->neighbor(1)->reconstructed()) {
        interior_face->neighbor(1)->reconstructed() = true;
        std::list<typename Triangulation::Vertex_handle> l1;
        get_boundary(interior_face->neighbor(1), interior_face->neighbor(1)->index(interior_face), l1);
        vertices.splice(vertices.end(), l1);
      } vertices.push_back(interior_face->vertex(2));
      if (interior_face->neighbor(0)->interior() && !interior_face->neighbor(0)->reconstructed()) {
        interior_face->neighbor(0)->reconstructed() = true;
        std::list<typename Triangulation::Vertex_handle> l0;
        get_boundary(interior_face->neighbor(0), interior_face->neighbor(0)->index(interior_face), l0);
        vertices.splice(vertices.end(), l0);
      } vertices.push_back(interior_face->vertex(1));

      // Find cutting vertices
      std::set<typename Triangulation::Vertex_handle> visited_vertices, repeated_vertices;
      for (auto const& current_vertex: vertices) {
        if (!visited_vertices.insert(current_vertex).second) repeated_vertices.insert(current_vertex);
      } visited_vertices.clear();

      // Cut and join rings in the correct order
      std::list<std::list<typename Triangulation::Vertex_handle>> rings;
      std::stack<std::list<typename Triangulation::Vertex_handle>> chains_stack;
      std::set<typename Triangulation::Vertex_handle> vertices_where_chains_begin;
      rings.push_back(std::list<typename Triangulation::Vertex_handle>());
      for (auto const& current_vertex: vertices) {
        
        // New chain
        if (repeated_vertices.count(current_vertex) > 0) {
          // Closed by itself
          if (rings.back().front() == current_vertex) {
            // Degenerate (insufficient vertices to be valid)
            if (rings.back().size() < 3) {
              rings.back().clear();
            } else {
              typename std::list<typename Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
              ++second_element;
              // Degenerate (zero area)
              if (rings.back().back() == *second_element) {
                rings.back().clear();
              }
              // Valid
              else {
                rings.push_back(std::list<typename Triangulation::Vertex_handle>());
              }
            }
          }
          // Open by itself
          else {
            // Closed with others in stack
            if (vertices_where_chains_begin.count(current_vertex)) {
              
              while (rings.back().front() != current_vertex) {
                rings.back().splice(rings.back().begin(), chains_stack.top());
                chains_stack.pop();
              } vertices_where_chains_begin.erase(current_vertex);
              // Degenerate (insufficient vertices to be valid)
              if (rings.back().size() < 3) {
                rings.back().clear();
              } else {
                typename std::list<typename Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
                ++second_element;
                // Degenerate (zero area)
                if (rings.back().back() == *second_element) {
                  rings.back().clear();
                }
                // Valid
                else {
                  rings.push_back(std::list<typename Triangulation::Vertex_handle>());
                }
              }
            }
            // Open
            else {
              // Not first chain
              if (repeated_vertices.count(rings.back().front()) > 0) {
                vertices_where_chains_begin.insert(rings.back().front());
              }
              chains_stack.push(std::list<typename Triangulation::Vertex_handle>());
              chains_stack.top().splice(chains_stack.top().begin(), rings.back());
            }
          }
        } rings.back().push_back(current_vertex);
      }
      // Final ring
      while (chains_stack.size() > 0) {
        rings.back().splice(rings.back().begin(), chains_stack.top());
        chains_stack.pop();
      }
      // Degenerate (insufficient vertices to be valid)
      if (rings.back().size() < 3) {
        rings.back().clear();
      } else {
        typename std::list<typename Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
        ++second_element;
        // Degenerate (zero area)
        if (rings.back().back() == *second_element) {
          rings.back().clear();
        }
      }
      
      // Remove last ring if too small (or empty)
      if (rings.back().size() < 3) {
        rings.pop_back();
      }
      
      // Start rings at the lexicographically smallest vertex
      for (auto& current_ring: rings) {
        typename std::list<typename Triangulation::Vertex_handle>::iterator smallest_vertex = current_ring.begin();
        for (typename std::list<typename Triangulation::Vertex_handle>::iterator current_vertex = current_ring.begin(); current_vertex != current_ring.end(); ++current_vertex) {
          if ((*current_vertex)->point() < (*smallest_vertex)->point()) smallest_vertex = current_vertex;
        } if (current_ring.back() != *smallest_vertex) {
          ++smallest_vertex;
          current_ring.splice(current_ring.begin(), current_ring, smallest_vertex, current_ring.end());
        }
      }

      // Make rings
      if (rings.size() == 0) continue;
      typename std::list<Polygon_2<Kernel, PolygonContainer>> rings_for_polygon;
      for (auto& current_ring: rings) {
        rings_for_polygon.push_back(Polygon_2<Kernel, PolygonContainer>());
        for (typename std::list<typename Triangulation::Vertex_handle>::reverse_iterator current_vertex = current_ring.rbegin(); current_vertex != current_ring.rend(); ++current_vertex) {
          rings_for_polygon.back().push_back(typename Kernel::Point_2((*current_vertex)->point().x(), (*current_vertex)->point().y()));
        }
      } typename std::list<Polygon_2<Kernel, PolygonContainer>>::iterator outer_ring;
      for (typename std::list<Polygon_2<Kernel, PolygonContainer>>::iterator current_ring = rings_for_polygon.begin(); current_ring != rings_for_polygon.end(); ++current_ring) {
        if (current_ring->orientation() == COUNTERCLOCKWISE) {
          outer_ring = current_ring;
          break;
        }
      } Polygon_with_holes_2<Kernel, PolygonContainer> new_polygon(*outer_ring);
      rings_for_polygon.erase(outer_ring);

      for (auto const& current_ring: rings_for_polygon) new_polygon.add_hole(current_ring);
      mp.add_polygon(new_polygon);
    }
  }

  // Erases the triangulation.
  void clear() {
    t.clear();
  }

  /// @}

  /// \name Access Functions
  /// @{

  Triangulation& triangulation() {
    return t;
  }

  Multipolygon_with_holes_2<Kernel, PolygonContainer> multipolygon() {
    return mp;
  }

  /// @}
  

protected:
  Triangulation t;
  Multipolygon_with_holes_2<Kernel, PolygonContainer> mp;
};

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H
