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
  Polygon_repair_2() {}

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
      face->added_to_list() = false;
    } 

    // Mark exterior as added to list and get interior triangles adjacent to it
    to_check_exterior.push_back(t.infinite_face());
    t.infinite_face()->added_to_list() = true;
    while (!to_check_exterior.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!to_check_exterior.front()->neighbor(neighbour)->added_to_list()) {
          if (!t.is_constrained(typename Triangulation::Edge(to_check_exterior.front(), neighbour))) {
            to_check_exterior.push_back(to_check_exterior.front()->neighbor(neighbour));
          } else {
            to_check_interior.push_back(to_check_exterior.front()->neighbor(neighbour));
          } to_check_exterior.front()->neighbor(neighbour)->added_to_list() = true;
        }
      } to_check_exterior.pop_front();
    }

    // Label region of front element of interior and exterior lists
    int current_polygon = 1, current_hole = -1;
    std::unordered_map<int, int> hole_nesting;
    while (!to_check_interior.empty() || !to_check_exterior.empty()) {

      // Interior
      if (!to_check_interior.empty()) {
        if (to_check_interior.front()->label() == 0) {
          to_check.push_back(to_check_interior.front());
          to_check_interior.front()->added_to_list() = true;
          while (!to_check.empty()) {
            to_check.front()->label() = current_polygon;
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (!t.is_constrained(typename Triangulation::Edge(to_check.front(), neighbour))) {
                to_check.push_back(to_check.front()->neighbor(neighbour));
              } else if (!to_check.front()->neighbor(neighbour)->added_to_list()) {
                to_check_exterior.push_back(to_check.front()->neighbor(neighbour));
                to_check_exterior_added_by.push_back(current_polygon);
              } to_check.front()->neighbor(neighbour)->added_to_list() = true;
            } to_check.pop_front();
          } ++current_polygon;
        } to_check_interior.pop_front();
      }

      // Exterior
      if (!to_check_exterior.empty()) {
        if (to_check_exterior.front()->label() == 0) {
          hole_nesting[current_hole] = to_check_exterior_added_by.front();
          to_check.push_back(to_check_exterior.front());
          to_check_exterior.front()->added_to_list() = true;
          while (!to_check.empty()) {
            to_check.front()->label() = current_hole;
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (!t.is_constrained(typename Triangulation::Edge(to_check.front(), neighbour))) {
                to_check.push_back(to_check.front()->neighbor(neighbour));
              } else if (!to_check.front()->neighbor(neighbour)->added_to_list()) {
                to_check_interior.push_back(to_check.front()->neighbor(neighbour));
              } to_check.front()->neighbor(neighbour)->added_to_list() = true;
            } to_check.pop_front();
          } --current_hole;
        } to_check_exterior.pop_front();
        to_check_exterior_added_by.pop_front();
      }

    }
  }

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
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
