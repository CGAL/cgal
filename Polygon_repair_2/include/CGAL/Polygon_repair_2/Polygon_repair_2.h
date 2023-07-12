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
  int number_of_polygons, number_of_holes;
  std::unordered_map<int, int> hole_nesting;
public:
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vertex_base;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
  typedef CGAL::Triangulation_face_base_with_repair_info_2<Kernel, Face_base> Face_base_with_repair_info;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info> Triangulation_data_structure;
  typedef CGAL::Exact_predicates_tag Tag; // assumed for now
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_Delaunay_triangulation;
  typedef Triangulation_with_odd_even_constraints_2<Constrained_Delaunay_triangulation> Triangulation;

  /// \name Creation
  Polygon_repair_2(): number_of_polygons(0), number_of_holes(0) {}

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

  // Label triangles in triangulation
  void label_triangulation() {
    std::list<typename Triangulation::Face_handle> to_check_exterior, to_check_interior, to_check;
    std::list<int> to_check_exterior_added_by;
    for (auto const face: t.all_face_handles()) {
      face->label() = 0;
      face->processed() = false;
    } 

    // Mark exterior as processed and put interior triangles adjacent to it in to_check_interior
    // these are used as starting points for the reconstruction
    // Note: exterior triangles are already labelled with 0
    to_check_exterior.push_back(t.infinite_face());
    t.infinite_face()->processed() = true; // processed means added to a list (to ensure elements are only added once)
    while (!to_check_exterior.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!to_check_exterior.front()->neighbor(neighbour)->processed()) {
          if (!t.is_constrained(typename Triangulation::Edge(to_check_exterior.front(), neighbour))) {
            to_check_exterior.push_back(to_check_exterior.front()->neighbor(neighbour));
          } else {
            to_check_interior.push_back(to_check_exterior.front()->neighbor(neighbour));
          } to_check_exterior.front()->neighbor(neighbour)->processed() = true;
        }
      } to_check_exterior.pop_front();
    }

    // Label region of front element of interior and exterior lists (alternating)
    while (!to_check_interior.empty() || !to_check_exterior.empty()) {

      // Interior triangle
      if (!to_check_interior.empty()) {
        if (to_check_interior.front()->label() == 0) { // label = 0 means not labelled yet
          to_check_interior.front()->label() = number_of_polygons+1;
          to_check.push_back(to_check_interior.front());
          while (!to_check.empty()) {
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (!t.is_constrained(typename Triangulation::Edge(to_check.front(), neighbour))) {
                if (to_check.front()->neighbor(neighbour)->label() == 0) {
                  to_check.front()->neighbor(neighbour)->label() = number_of_polygons+1;
                  to_check.push_back(to_check.front()->neighbor(neighbour));
                  to_check.front()->neighbor(neighbour)->processed() = true;
                } 
              } else { // constrained -> exterior or hole
                if (!to_check.front()->neighbor(neighbour)->processed()) { // only add holes (not processed) to list once
                  to_check_exterior.push_back(to_check.front()->neighbor(neighbour));
                  to_check.front()->neighbor(neighbour)->processed() = true;
                  to_check_exterior_added_by.push_back(number_of_polygons+1);
                } 
              }
            } to_check.pop_front();
          } ++number_of_polygons;
        } to_check_interior.pop_front();
      }

      // Exterior triangle (hole)
      if (!to_check_exterior.empty()) {
        if (to_check_exterior.front()->label() == 0) { // label = 0 means not labelled yet
          to_check_exterior.front()->label() = -number_of_holes+1;
          to_check.push_back(to_check_exterior.front());
          hole_nesting[-number_of_holes+1] = to_check_exterior_added_by.front(); // record nesting of current hole
          while (!to_check.empty()) {
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (!t.is_constrained(typename Triangulation::Edge(to_check.front(), neighbour))) {
                if (to_check.front()->neighbor(neighbour)->label() == 0) {
                  to_check.front()->neighbor(neighbour)->label() = number_of_polygons+1;
                  to_check.push_back(to_check.front()->neighbor(neighbour));
                  to_check.front()->neighbor(neighbour)->processed() = true;
                }
              } else { // constrained -> interior
                if (!to_check.front()->neighbor(neighbour)->processed()) { // interior triangles only added once
                  to_check_interior.push_back(to_check.front()->neighbor(neighbour));
                  to_check.front()->neighbor(neighbour)->processed() = true;
                }
              }
            } to_check.pop_front();
          } ++number_of_holes;
        } to_check_exterior.pop_front();
        to_check_exterior_added_by.pop_front();
      }

    } std::cout << number_of_polygons << " polygons with " << number_of_holes << " holes in triangulation" << std::endl;
  }

  // Reconstruct ring boundary starting from an edge (face + opposite vertex) that is part of it
  void reconstruct_ring(Polygon_2<Kernel, PolygonContainer>& ring,
                        typename Triangulation::Face_handle face_adjacent_to_boundary, int opposite_vertex) {
    typename Triangulation::Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    do {
      typename Triangulation::Vertex_handle pivot_vertex = current_face->vertex(current_face->cw(current_opposite_vertex));
      // std::cout << "Adding point " << pivot_vertex->point() << std::endl;
      ring.push_back(pivot_vertex->point());
      typename Triangulation::Face_circulator fc = t.incident_faces(pivot_vertex, current_face);
      do {
        ++fc;
      } while (fc->label() != current_face->label());
      current_face = fc;
      current_opposite_vertex = fc->cw(fc->index(pivot_vertex));
    } while (current_face != face_adjacent_to_boundary ||
             current_opposite_vertex != opposite_vertex);
  }

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
    std::vector<Polygon_2<Kernel, PolygonContainer>> polygons, holes;
    polygons.reserve(number_of_polygons);
    holes.reserve(number_of_holes);

    for (auto const face: t.all_face_handles()) {
      face->processed() = false;
    } for (auto const &face: t.finite_face_handles()) {
      if (face->label() == 0) continue; // exterior triangle
      if (face->processed()) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
        if (face->label() != face->neighbor(opposite_vertex)->label()) {

          if (face->label() > 0) {
            polygons.emplace_back();
            reconstruct_ring(polygons.back(), face, opposite_vertex);
          } else {
            holes.emplace_back();
            reconstruct_ring(holes.back(), face, opposite_vertex);
          } 

          std::list<typename Triangulation::Face_handle> to_check;
          to_check.push_back(face);
          while (!to_check.empty()) {
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (to_check.front()->label() == to_check.front()->neighbor(neighbour)->label() &&
                  !to_check.front()->neighbor(neighbour)->processed()) {
                to_check.push_back(to_check.front()->neighbor(neighbour));
              } 
            } to_check.front()->processed() = true;
            to_check.pop_front();
          } 

          break;
        }
      }
    }

    for (auto& polygon: polygons) {
      // std::cout << "Adding polygon " << polygon << std::endl;
      mp.add_polygon(Polygon_with_holes_2<Kernel, PolygonContainer>(polygon));
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
