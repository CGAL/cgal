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

#ifndef CGAL_POLYGON_REPAIR_H
#define CGAL_POLYGON_REPAIR_H

#include <CGAL/license/Polygon_repair.h>

#include <list>
#include <vector>
#include <set>
#include <unordered_set>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_repair/Even_odd_rule.h>

#include <CGAL/Polygon_repair/internal/Triangulation_face_base_with_repair_info_2.h>
#include <CGAL/Polygon_repair/internal/Triangulation_with_even_odd_constraints_2.h>

namespace CGAL {

namespace Polygon_repair {

#ifndef DOXYGEN_RUNNING
template <class Kernel, class Container>
class Polygon_repair;
#endif

/// \ingroup PkgPolygonRepairFunctions
/// repairs polygon `p` using the given rule
/// \tparam Kernel parameter of the input and output polygons
/// \tparam Container parameter of the input and output polygons
///  \tparam Rule must be `Even_odd_rule`
template <class Kernel, class Container, class Rule = Even_odd_rule>
Multipolygon_with_holes_2<Kernel, Container> repair(const Polygon_2<Kernel, Container>& p , Rule = Rule())
{
  static_assert(std::is_same_v<Rule,Even_odd_rule>);
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container> pr;
  pr.add_to_triangulation_even_odd(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_even_odd();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

/// \ingroup PkgPolygonRepairFunctions
/// repairs polygon with holes `p` using the given rule
/// \tparam Kernel parameter of the input and output polygons
/// \tparam Container parameter of the input and output polygons
///  \tparam Rule must be `Even_odd_rule`
template <class Kernel, class Container, class Rule = Even_odd_rule>
Multipolygon_with_holes_2<Kernel, Container> repair(const Polygon_with_holes_2<Kernel, Container>& p, Rule = Rule())
{
  static_assert(std::is_same_v<Rule,Even_odd_rule>);
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container> pr;
  pr.add_to_triangulation_even_odd(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_even_odd();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

/// \ingroup PkgPolygonRepairFunctions
/// repairs multipolygon with holes `p` using the given rule
/// \tparam Kernel parameter of the input and output polygons
/// \tparam Container parameter of the input and output polygons
///  \tparam Rule must be `Even_odd_rule`
template <class Kernel, class Container, class Rule = Even_odd_rule>
Multipolygon_with_holes_2<Kernel, Container> repair(const Multipolygon_with_holes_2<Kernel, Container>& p, Rule = Rule())
{
  static_assert(std::is_same_v<Rule,Even_odd_rule>);
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container> pr;
  pr.add_to_triangulation_even_odd(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_even_odd();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

template <class Kernel, class Container>
bool is_valid(const Polygon_2<Kernel, Container>& polygon) {
  if (polygon.vertices().size() < 3) {
    std::cout << "Invalid: less than 3 vertices" << std::endl;
    return false;
  } for (auto const& edge: polygon.edges()) {
    if (edge.source() == edge.target()) {
      std::cout << "Invalid: duplicate vertices" << std::endl;
      return false;
    }
  } if (!polygon.is_simple()) {
    std::cout << "Invalid: not simple" << std::endl;
    return false;
  } return true;
}

template <class Kernel, class Container>
bool is_valid(const Polygon_with_holes_2<Kernel, Container>& polygon) {

  // Validate outer boundary
  for (auto const& edge: polygon.outer_boundary().edges()) {
    if (edge.source() == edge.target()) {
      std::cout << "Invalid: duplicate vertices in outer boundary" << std::endl;
      return false;
    }
  } if (!polygon.outer_boundary().is_simple()) {
    std::cout << "Invalid: outer boundary not simple" << std::endl;
    return false;
  }

  // Validate holes
  for (auto const& hole: polygon.holes()) {
    for (auto const& edge: hole.edges()) {
      if (edge.source() == edge.target()) {
        std::cout << "Invalid: duplicate vertices in hole" << std::endl;
        return false;
      }
    } if (!hole.is_simple()) {
      std::cout << "Invalid: hole not simple" << std::endl;
      return false;
    }
  }

  // Create triangulation of outer boundary
  typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation vt;
  for (auto const& edge: polygon.outer_boundary().edges()) {
    try {
      vt.insert_constraint(edge.source(), edge.target());
    } catch (typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Intersection_of_constraints_exception ice) {
      std::cout << "Invalid: intersection in outer boundary" << std::endl;
      return false;
    }
  } if (vt.number_of_faces() == 0) {
    std::cout << "Invalid: no outer boundary" << std::endl;
    return false;
  } for (auto const face: vt.all_face_handles()) {
    face->label() = 0;
    face->processed() = false;
  } std::list<typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Face_handle> to_check;
  std::list<int> to_check_added_by;
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by); // exterior
  int regions = 0, holes = 0;
  while (!to_check.empty()) {
    if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
      if (to_check_added_by.front() < 0) {
        CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by);
        ++regions;
      } else {
        CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by);
        ++holes;
      }
    } to_check.pop_front();
    to_check_added_by.pop_front();
  } CGAL_assertion(regions == 1 && holes == 0);

  // Hole nesting
  for (auto const& hole: polygon.holes()) {
    for (auto const& vertex: hole.vertices()) {
      typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Locate_type lt;
      int li;
      typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Face_handle f = vt.locate(vertex, lt, li);
      if (lt == CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Locate_type::FACE && f->label() != 1) {
        std::cout << "Invalid: hole (partly) outside outer boundary" << std::endl;
        return false;
      }
    }
    for (auto const& edge: hole.edges()) {
      try {
        vt.insert_constraint(edge.source(), edge.target());
      } catch (typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Intersection_of_constraints_exception ice) {
        std::cout << "Invalid: hole (partly) outside outer boundary" << std::endl;
        return false;
      }
    }
  }

  // Connected interior
  for (auto const face: vt.all_face_handles()) {
    face->label() = 0;
    face->processed() = false;
  } to_check.clear();
  to_check_added_by.clear();
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by); // exterior
  regions = 0;
  holes = 0;
  while (!to_check.empty()) {
    if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
      if (to_check_added_by.front() < 0) {
        CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by);
        ++regions;
      } else {
        CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by);
        ++holes;
      }
    } to_check.pop_front();
    to_check_added_by.pop_front();
  } if (regions != 1) {
    std::cout << "Invalid: disconnected interior" << std::endl;
    return false;
  } CGAL_assertion(holes == polygon.number_of_holes());

  return true;
}

template <class Kernel, class Container>
bool is_valid(const Multipolygon_with_holes_2<Kernel, Container>& multipolygon) {

  // Validate polygons
  for (auto const& polygon: multipolygon.polygons_with_holes()) {
    if (!is_valid(polygon)) return false;
  }

  typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation vt;
  typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Locate_type lt;
  int li;
  for (auto const& polygon: multipolygon.polygons_with_holes()) {


    if (vt.number_of_faces() > 0) {

      // Relabel
      for (auto const face: vt.all_face_handles()) {
        face->label() = 0;
        face->processed() = false;
      } std::list<typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Face_handle> to_check;
      std::list<int> to_check_added_by;
      CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by); // exterior
      int regions = 0, holes = 0;
      while (!to_check.empty()) {
        if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
          if (to_check_added_by.front() < 0) {
            CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by);
            ++regions;
          } else {
            CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by);
            ++holes;
          }
        } to_check.pop_front();
        to_check_added_by.pop_front();
      }

      // Test vertices in labeled triangulation
      for (auto const& vertex: polygon.outer_boundary().vertices()) {
        typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Face_handle f = vt.locate(vertex, lt, li);
        if (lt == CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Locate_type::FACE && f->label() > 0) {
          std::cout << "Invalid: (partly) overlapping polygons" << std::endl;
          return false;
        }
      }
      for (auto const& hole: polygon.holes()) {
        for (auto const& vertex: hole.vertices()) {
          typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Face_handle f = vt.locate(vertex, lt, li);
          if (lt == CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Locate_type::FACE && f->label() > 0) {
            std::cout << "Invalid: (partly) overlapping polygons" << std::endl;
            return false;
          }
        }
      }

    }

    // Insert constraints while checking for intersections
    for (auto const& edge: polygon.outer_boundary().edges()) {
      try {
        vt.insert_constraint(edge.source(), edge.target());
      } catch (typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Intersection_of_constraints_exception ice) {
        std::cout << "Invalid: (partly) overlapping polygons" << std::endl;
        return false;
      }
    }
    for (auto const& hole: polygon.holes()) {
      for (auto const& edge: hole.edges()) {
        try {
          vt.insert_constraint(edge.source(), edge.target());
        } catch (typename CGAL::Polygon_repair::Polygon_repair<Kernel, Container>::Validation_triangulation::Intersection_of_constraints_exception ice) {
          std::cout << "Invalid: (partly) overlapping polygons" << std::endl;
          return false;
        }
      }
    }
  }

  return true;
}

#ifndef DOXYGEN_RUNNING

template <class Kernel, class Container = std::vector<typename Kernel::Point_2>>
class Polygon_repair {
public:
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Vertex_base = CGAL::Triangulation_vertex_base_2<Kernel>;
  using Face_base = CGAL::Constrained_triangulation_face_base_2<Kernel>;
  using Face_base_with_repair_info = internal::Triangulation_face_base_with_repair_info_2<Kernel, Face_base>;
  using Triangulation_data_structure = CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info>;
  using Tag = typename std::conditional<std::is_floating_point<FT>::value,
                                        CGAL::Exact_predicates_tag,
                                        CGAL::Exact_intersections_tag>::type;
  using Constrained_Delaunay_triangulation = CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag>;
  using Triangulation = internal::Triangulation_with_even_odd_constraints_2<Constrained_Delaunay_triangulation>;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Face_handle = typename Triangulation::Face_handle;
  using Face_circulator = typename Triangulation::Face_circulator;
  using Edge = typename Triangulation::Edge;

  using Edge_map = typename std::conditional<std::is_floating_point<FT>::value,
                                             std::unordered_set<std::pair<Point_2, Point_2>,
                                                                boost::hash<std::pair<Point_2, Point_2>>>,
                                             std::set<std::pair<Point_2, Point_2>>>::type;
  using Vertex_map = typename std::conditional<std::is_floating_point<FT>::value,
                                               std::unordered_map<Point_2, Vertex_handle>,
                                               std::map<Point_2, Vertex_handle>>::type;

  using Validation_tag = CGAL::No_constraint_intersection_tag;
  using Validation_triangulation = CGAL::Constrained_triangulation_2<Kernel, Triangulation_data_structure, Validation_tag>;

  using Polygon_2 = CGAL::Polygon_2<Kernel, Container>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel, Container>;
  using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel, Container>;

  struct Polygon_less {
    bool operator()(const Polygon_2& pa, const Polygon_2& pb) const {
      typename Polygon_2::Vertex_iterator va = pa.vertices_begin();
      typename Polygon_2::Vertex_iterator vb = pb.vertices_begin();
      while (va != pa.vertices_end() && vb != pb.vertices_end()) {
        if (*va != *vb) return *va < *vb;
        ++va;
        ++vb;
      }
      if (vb == pb.vertices_end()) return false;
      return true;
    }
  };

  struct Polygon_with_holes_less {
    Polygon_less pl;
    bool operator()(const Polygon_with_holes_2& pa, const Polygon_with_holes_2& pb) const {
      if (pl(pa.outer_boundary(), pb.outer_boundary())) return true;
      if (pl(pb.outer_boundary(), pa.outer_boundary())) return false;
      typename Polygon_with_holes_2::Hole_const_iterator ha = pa.holes_begin();
      typename Polygon_with_holes_2::Hole_const_iterator hb = pb.holes_begin();
      while (ha != pa.holes_end() && hb != pb.holes_end()) {
        if (pl(*ha, *hb)) return true;
        if (pl(*hb, *ha)) return false;
      }
      if (hb == pb.holes_end()) return false;
      return true;
    }
  };

  /// \name Creation
  Polygon_repair() : number_of_polygons(0), number_of_holes(0) {}

  /// \name Modifiers
  /// @{

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_even_odd(const Polygon_2& polygon) {

    // Get unique edges
    for (auto const& edge: polygon.edges()) {
      if (edge.source() == edge.target()) continue;
      std::pair<Point_2, Point_2> pair = (edge.source() < edge.target())?
      std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
      auto inserted = unique_edges.insert(pair);
      if (!inserted.second) unique_edges.erase(inserted.first);
    }

    // Insert vertices
    Vertex_map vertices;
    std::vector<std::pair<Vertex_handle, Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      Vertex_handle first_vertex, second_vertex;
      typename Vertex_map::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.even_odd_insert_constraint(edge.first, edge.second);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_even_odd(const Polygon_with_holes_2& polygon) {

    // Get unique edges
    for (auto const& edge: polygon.outer_boundary().edges()) {
      if (edge.source() == edge.target()) continue;
      std::pair<Point_2, Point_2> pair = (edge.source() < edge.target())?
      std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
      auto inserted = unique_edges.insert(pair);
      if (!inserted.second) unique_edges.erase(inserted.first);
    }
    for (auto const& hole: polygon.holes()) {
      for (auto const& edge: hole.edges()) {
        if (edge.source() == edge.target()) continue;
        std::pair<Point_2, Point_2> pair = (edge.source() < edge.target())?
        std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
        auto inserted = unique_edges.insert(pair);
        if (!inserted.second) unique_edges.erase(inserted.first);
      }
    }

    // Insert vertices
    Vertex_map vertices;
    std::vector<std::pair<Vertex_handle, Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      Vertex_handle first_vertex, second_vertex;
      typename Vertex_map::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.even_odd_insert_constraint(edge.first, edge.second);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_even_odd(const Multipolygon_with_holes_2& multipolygon) {

    // Get unique edges
    for (auto const& polygon: multipolygon.polygons_with_holes()) {
      for (auto const& edge: polygon.outer_boundary().edges()) {
        if (edge.source() == edge.target()) continue;
        std::pair<Point_2, Point_2> pair = (edge.source() < edge.target())?
        std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
        auto inserted = unique_edges.insert(pair);
        if (!inserted.second) unique_edges.erase(inserted.first);
      }
      for (auto const& hole: polygon.holes()) {
        for (auto const& edge: hole.edges()) {
          if (edge.source() == edge.target()) continue;
          std::pair<Point_2, Point_2> pair = (edge.source() < edge.target())?
          std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
          auto inserted = unique_edges.insert(pair);
          if (!inserted.second) unique_edges.erase(inserted.first);
        }
      }
    }

    // Insert vertices
    Vertex_map vertices;
    std::vector<std::pair<Vertex_handle, Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      Vertex_handle first_vertex, second_vertex;
      typename Vertex_map::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.even_odd_insert_constraint(edge.first, edge.second);
    }
  }

  // Label a region of adjacent triangles without passing through constraints
  // adjacent triangles that involve passing through constraints are added to to_check
  template <class T>
  static void label_region(T& tt, Face_handle face, int label,
                           std::list<Face_handle>& to_check,
                           std::list<int>& to_check_added_by) {
    // std::cout << "Labelling region with " << label << std::endl;
    std::list<Face_handle> to_check_in_region;
    face->label() = label;
    to_check_in_region.push_back(face);
    face->processed() = true; // processed means added to a list (to ensure elements are only added once)

    while (!to_check_in_region.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!tt.is_constrained(Edge(to_check_in_region.front(), neighbour))) {
          if (to_check_in_region.front()->neighbor(neighbour)->label() == 0) { // unlabeled
            to_check_in_region.front()->neighbor(neighbour)->label() = label;
            to_check_in_region.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        } else { // constrained
          if (!to_check_in_region.front()->neighbor(neighbour)->processed()) { // not added to to_check
            to_check.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_added_by.push_back(label);
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        }
      } to_check_in_region.pop_front();
    }
  }

  // Label triangles in triangulation
  void label_triangulation_even_odd() {

    // Simplify collinear edges (gets rid of order dependency)
    for (auto vertex: t.all_vertex_handles()) {
      typename Triangulation::Edge_circulator first_edge = t.incident_edges(vertex);
      typename Triangulation::Edge_circulator current_edge = first_edge;
      std::vector<Edge> incident_constrained_edges;
      do {
        if (t.is_constrained(*current_edge)) {
          incident_constrained_edges.push_back(*current_edge);
        } ++current_edge;
      } while (current_edge != first_edge);
      if (incident_constrained_edges.size() == 2) {
        Point_2 v1 = incident_constrained_edges.front().first->vertex(incident_constrained_edges.front().first->ccw(incident_constrained_edges.front().second))->point();
        Point_2 v2 = incident_constrained_edges.back().first->vertex(incident_constrained_edges.back().first->ccw(incident_constrained_edges.back().second))->point();
        if (CGAL::collinear(v1, vertex->point(), v2)) {
          // std::cout << "Collinear points" << std::endl;
          // std::cout << "v1: " << v1 << std::endl;
          // std::cout << "in: " << vertex->point() << std::endl;
          // std::cout << "v2: " << v2 << std::endl;
          t.remove_incident_constraints(vertex);
          t.remove(vertex);
          t.insert_constraint(v1, v2);
        }
      }
    }

    // Init labels
    for (auto const face: t.all_face_handles()) {
      face->label() = 0;
      face->processed() = false;
    }

    // Label exterior with label -1, marking it as processed and
    // putting interior triangles adjacent to it in to_check
    std::list<Face_handle> to_check;
    std::list<int> to_check_added_by;
    label_region(t, t.infinite_face(), -1, to_check, to_check_added_by);

    // Label region of front element to_check list
    while (!to_check.empty()) {

      if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
        if (to_check_added_by.front() < 0) {
          label_region(t, to_check.front(), number_of_polygons+1, to_check, to_check_added_by);
          ++number_of_polygons;
        } else {
          label_region(t, to_check.front(), -(number_of_holes+2), to_check, to_check_added_by);
          ++number_of_holes;
        }
      } to_check.pop_front();
      to_check_added_by.pop_front();

    } // std::cout << number_of_polygons << " polygons with " << number_of_holes << " holes in triangulation" << std::endl;
  }

  // Reconstruct ring boundary starting from an edge (face + opposite vertex) that is part of it
  void reconstruct_ring(std::list<Point_2>& ring,
                        Face_handle face_adjacent_to_boundary,
                        int opposite_vertex) {
    // std::cout << "Reconstructing ring for face " << face_adjacent_to_boundary->label() << "..." << std::endl;

    // Create ring
    Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    do {
      CGAL_assertion(current_face->label() == face_adjacent_to_boundary->label());
      current_face->processed() = true;
      Vertex_handle pivot_vertex = current_face->vertex(current_face->cw(current_opposite_vertex));
      // std::cout << "\tAdding point " << pivot_vertex->point() << std::endl;
      ring.push_back(pivot_vertex->point());
      Face_circulator fc = t.incident_faces(pivot_vertex, current_face);
      do {
        ++fc;
      } while (fc->label() != current_face->label());
      current_face = fc;
      current_opposite_vertex = fc->cw(fc->index(pivot_vertex));
    } while (current_face != face_adjacent_to_boundary ||
             current_opposite_vertex != opposite_vertex);

    // Start at lexicographically smallest vertex
    typename std::list<Point_2>::iterator smallest_vertex = ring.begin();
    for (typename std::list<Point_2>::iterator current_vertex = ring.begin();
         current_vertex != ring.end(); ++current_vertex) {
      if (*current_vertex < *smallest_vertex) smallest_vertex = current_vertex;
    }
    if (ring.front() != *smallest_vertex) {
      ring.splice(ring.begin(), ring, smallest_vertex, ring.end());
    }
  }

  // Reconstruct multipolygon based on the triangles labeled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
    std::vector<Polygon_2> polygons; // outer boundaries
    std::vector<std::set<Polygon_2, Polygon_less>> holes; // holes are ordered (per polygon)
    polygons.resize(number_of_polygons);
    holes.resize(number_of_polygons);

    for (auto const face: t.all_face_handles()) {
      face->processed() = false;
    }
    for (auto const &face: t.finite_face_handles()) {
      if (face->label() < 1) continue; // exterior triangle
      if (face->processed()) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
        if (face->label() == face->neighbor(opposite_vertex)->label()) continue; // not adjacent to boundary

        // Reconstruct ring
        std::list<Point_2> ring;
        reconstruct_ring(ring, face, opposite_vertex);

        // Put ring in polygons
        Polygon_2 polygon;
        polygon.reserve(ring.size());
        polygon.insert(polygon.vertices_end(),ring.begin(), ring.end());
        // std::cout << "Reconstructed ring for polygon " << face->label() << " with ccw? " << (polygon.orientation() == CGAL::COUNTERCLOCKWISE) << std::endl;
        if (polygon.orientation() == CGAL::COUNTERCLOCKWISE) {
          polygons[face->label()-1] = std::move(polygon);
        } else {
          holes[face->label()-1].insert(std::move(polygon));
        } break;
      }
    }

    // Create polygons with holes and put in multipolygon
    std::set<Polygon_with_holes_2, Polygon_with_holes_less> ordered_polygons;
    for (std::size_t i = 0; i < polygons.size(); ++i) {
      ordered_polygons.insert(Polygon_with_holes_2(std::move(polygons[i]), std::make_move_iterator(holes[i].begin()), std::make_move_iterator(holes[i].end())));
    }
    for (auto const& polygon: ordered_polygons) {
      // std::cout << "Adding polygon " << polygon << std::endl;
      mp.add_polygon_with_holes(std::move(polygon));
    }
  }

  // Erases the triangulation.
  void clear() {
    t.clear();
    unique_edges.clear();
    mp.clear();
    number_of_polygons = 0;
    number_of_holes = 0;
    search_start = Triangulation::Face_handle();
  }

  /// @}

  /// \name Access Functions
  /// @{

  Triangulation& triangulation() {
    return t;
  }

  const Multipolygon_with_holes_2& multipolygon() {
    return mp;
  }

  /// @}


protected:
  Triangulation t;
  Edge_map unique_edges;
  Multipolygon_with_holes_2 mp;
  int number_of_polygons, number_of_holes;
  Face_handle search_start;
};

#endif // DOXYGEN_RUNNING

} // namespace Polygon_repair
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_H
