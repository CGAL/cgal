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
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polygon_repair/Even_odd_rule.h>
#include <CGAL/Polygon_repair/Non_zero_rule.h>
#include <CGAL/Polygon_repair/Union_rule.h>
#include <CGAL/Polygon_repair/Intersection_rule.h>

#include <CGAL/Polygon_repair/internal/Triangulation_face_base_with_repair_info_2.h>
#include <CGAL/Polygon_repair/internal/Triangulation_with_even_odd_constraints_2.h>
#include <CGAL/Polygon_repair/Winding.h>
#include <CGAL/Polygon_repair/Boolean.h>

namespace CGAL {

namespace Polygon_repair {

#ifndef DOXYGEN_RUNNING
template <class Kernel, class Container>
class Polygon_repair;
#endif

/// \ingroup PkgPolygonRepairFunctions
/// repairs polygon `p` using the given rule
/// \tparam Kernel parameter of the input and output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
///  \tparam Rule must be `Even_odd_rule` or `Non_zero_rule`
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
/// \tparam Kernel parameter of the input and output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
///  \tparam Rule must be `Even_odd_rule` or `Non_zero_rule`
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


template <class Kernel, class Container>
Multipolygon_with_holes_2<Kernel, Container> repair(const Polygon_2<Kernel, Container>& p, Non_zero_rule)
{
  Winding<Kernel> winding;
  winding.insert(p);
  winding.label();
  winding.label_domains();
  return winding();
}


template <class Kernel, class Container>
Multipolygon_with_holes_2<Kernel, Container> repair(const Polygon_with_holes_2<Kernel, Container>& p, Non_zero_rule)
{
  Winding<Kernel> winding;
  winding.insert(p);
  winding.label();
  winding.label_domains();
  return winding();
}


/// \ingroup PkgPolygonRepairFunctions
/// repairs multipolygon with holes `p` using the given rule
/// \tparam Kernel parameter of the input and output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
/// \tparam Rule may be any \ref PkgPolygonRepairRules
/// \pre If the rule is the `Union_rule` or `Non_zero_rule`, each polygon with hole must be free of self-intersections,
///      the outer boundary of each polygon with holes must be counterclockwise  and the  holes clockwise oriented.
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
Multipolygon_with_holes_2<Kernel, Container> repair(const Multipolygon_with_holes_2<Kernel, Container>& p, Non_zero_rule)
{
  Winding<Kernel> winding;
  winding.insert(p);
  winding.label();
  winding.label_domains();
  return winding();
}


template <class Kernel, class Container>
Multipolygon_with_holes_2<Kernel, Container> repair(const Multipolygon_with_holes_2<Kernel, Container>& p, Union_rule)
{
  struct Larger_than_zero {
    bool operator()(int i) const
    {
      return i > 0;
    }
  };

  CGAL::Polygon_repair::Boolean<Kernel,Container> bops;
  bops.insert(p);
  bops.mark_domains();
  Larger_than_zero ltz;
  return bops(ltz);
}


template <class Kernel, class Container>
Multipolygon_with_holes_2<Kernel, Container> repair(const Multipolygon_with_holes_2<Kernel, Container>& p, Intersection_rule)
{
  struct Equal  {
    int val;
    Equal(int val)
    : val(val)
    {}

    bool operator()(int i) const
    {
      return i == val;
    }
  };

 CGAL::Polygon_repair::Boolean<Kernel,Container> bops;
  bops.insert(p);
  bops.mark_domains();
  Equal equal(p.number_of_polygons_with_holes());
  return bops(equal);
}

/*
template <class Kernel, class Container>
Multipolygon_with_holes_2<Kernel, Container> repair(const Multipolygon_with_holes_2<Kernel, Container>& p, Non_zero_rule rule)
{
  static_assert(std::is_same_v<Rule,Even_odd_rule>);
  CGAL::Polygon_repair::Polygon_repair<Kernel, Container> pr;
  pr.add_to_triangulation_even_odd(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_even_odd();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}
*/

template <class Kernel, class Container>
bool is_valid(const Polygon_2<Kernel, Container>& polygon) {
  if (polygon.vertices().size() < 3) {
    std::cerr << "Invalid: fewer than 3 vertices" << std::endl;
    return false;
  }
  for (auto const& edge: polygon.edges()) {
    if (edge.source() == edge.target()) {
      std::cerr << "Invalid: degenerate edge" << std::endl;
      return false;
    }
  }
  if (!polygon.is_simple()) {
    std::cerr << "Invalid: not simple" << std::endl;
    return false;
  }
  return true;
}

template <class Kernel, class Container>
bool is_valid(const Polygon_with_holes_2<Kernel, Container>& polygon)
{
  using PR = CGAL::Polygon_repair::Polygon_repair<Kernel, Container>;
  using VTr = typename PR::Validation_triangulation;
  using Vertex_handle = typename VTr::Vertex_handle;
  using Edge = typename VTr::Edge;
  using Face_handle = typename VTr::Face_handle;
  using Context = typename VTr::Context;
  using Constraint_id = typename VTr::Constraint_id;

  // Validate outer boundary
  if (!is_valid(polygon.outer_boundary())) {
    std::cerr << "Invalid: invalid outer hull" << std::endl;
    return false;
  }

  if (polygon.outer_boundary().orientation() != CGAL::COUNTERCLOCKWISE) {
    std::cerr << "Invalid: wrong outer hull orientation" << std::endl;
    return false;
  }

  // Validate holes
  for (auto const& hole: polygon.holes()) {
    if (!is_valid(hole)) {
      std::cerr << "Invalid: invalid hole" << std::endl;
      return false;
    }
    if (hole.orientation() != CGAL::CLOCKWISE) {
      std::cerr << "Invalid: wrong hole orientation" << std::endl;
      return false;
    }
  }

  // Create constrained triangulation of outer boundary
  std::set<Constraint_id> outer_hull_constraints;

  VTr vt;
  for (auto const& edge: polygon.outer_boundary().edges()) {
    try {
      Constraint_id cid = vt.insert_constraint(edge.source(), edge.target());
      outer_hull_constraints.insert(cid);
    } catch (typename VTr::Intersection_of_constraints_exception&) {
      std::cerr << "Invalid: intersection in outer boundary" << std::endl;
      return false;
    }
  }

  CGAL_postcondition(outer_hull_constraints.size() == polygon.outer_boundary().edges().size());

  if (vt.number_of_faces() == 0) {
    std::cerr << "Invalid: no outer boundary" << std::endl;
    return false;
  }

  // Add all hole boundaries as constraints
  std::set<Constraint_id> hole_constraints;

  for (auto const& hole: polygon.holes()) {
    for (auto const& edge: hole.edges()) {
      try {
        Constraint_id cid = vt.insert_constraint(edge.source(), edge.target());
        hole_constraints.insert(cid);
      } catch (typename VTr::Intersection_of_constraints_exception&) {
        std::cerr << "Invalid: hole boundary intersects something" << std::endl;
        return false;
      }
    }
  }

  for (auto const face: vt.all_face_handles()) {
    face->label() = 0;
    face->processed() = false;
  }

  // A hole border is valid only if its two incident faces are inside the outer hull.
  // We want to cross through hole constraint, but not outer hull constraint.
  struct Outer_hull_constraint_check
  {
    Outer_hull_constraint_check(const std::set<Constraint_id>& outer_hull_constraints)
      : m_outer_hull_constraints(outer_hull_constraints)
    { }

    bool operator()(const Edge& e, const VTr& vt) const
    {
      if (!vt.is_constrained(e)) {
        return false;
      }
      Vertex_handle va = e.first->vertex((e.second + 1)%3);
      Vertex_handle vb = e.first->vertex((e.second + 2)%3);
      for (const Context& c : vt.contexts(va, vb)) {
        if (m_outer_hull_constraints.count(c.id()) != 0) {
          return true;
        }
      }
      return false;
    };

    std::set<Constraint_id> m_outer_hull_constraints;
  };

  Outer_hull_constraint_check constraint_checker(outer_hull_constraints);

  std::list<typename VTr::Face_handle> to_check;
  std::list<int> to_check_added_by;
  PR::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by, constraint_checker); // exterior
  int regions = 0, holes = 0;
  while (!to_check.empty()) {
    if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
      if (to_check_added_by.front() < 0) {
        PR::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by, constraint_checker);
        ++regions;
      } else {
        PR::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by, constraint_checker);
        ++holes;
      }
    }
    to_check.pop_front();
    to_check_added_by.pop_front();
  }

  CGAL_assertion(regions == 1 && holes == 0);

  // Now, look at the hole edge constraints
  for (const Edge& e : vt.all_edges()) {
    if (!vt.is_constrained(e)) {
      continue;
    }

    Vertex_handle va = e.first->vertex((e.second + 1)%3);
    Vertex_handle vb = e.first->vertex((e.second + 2)%3);
    unsigned int hits = 0;
    for (const Context& c : vt.contexts(va, vb)) {
      if (hole_constraints.count(c.id()) != 0) {
        ++hits;
      }
    }

    if (hits == 0) { // edge is not a hole constraint
      continue;
    }

    if (hits >= 2) { // edge appears in multiple holes
      CGAL_assertion(polygon.number_of_holes() > 1); // holes are known to be simple here
      std::cerr << "Invalid: hole edge appears multiple times" << std::endl;
      return false;
    }

    // from now on, the edge is hole boundary appearing only once

    Face_handle fh = e.first;
    Face_handle nfh = e.first->neighbor(e.second);

    // hole boundary should be entirely within the area delimited by the polygon's outer boundary
    if (fh->label() != 1 || nfh->label() != 1) {
      std::cerr << "Invalid: outward hole edge" << std::endl;
      return false;
    }
  }

  // Connected interior detection: flood again, but take into account hole constraints too.
  // There should be a single region
  for (auto const face: vt.all_face_handles()) {
    face->label() = 0;
    face->processed() = false;
  }
  to_check.clear();
  to_check_added_by.clear();

  PR::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by); // exterior

  regions = 0;
  holes = 0;
  while (!to_check.empty()) {
    if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
      if (to_check_added_by.front() < 0) {
        PR::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by);
        ++regions;
      } else {
        PR::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by);
        ++holes;
      }
    }
    to_check.pop_front();
    to_check_added_by.pop_front();
  }

  if (regions != 1) {
    std::cerr << "Invalid: disconnected interior" << std::endl;
    return false;
  }

  CGAL_assertion(holes == static_cast<int>(polygon.number_of_holes()));

  return true;
}

template <class Kernel, class Container>
bool is_valid(const Multipolygon_with_holes_2<Kernel, Container>& multipolygon)
{
  using PR = CGAL::Polygon_repair::Polygon_repair<Kernel, Container>;
  using VTr = typename PR::Validation_triangulation;
  using Vertex_handle = typename VTr::Vertex_handle;
  using Edge = typename VTr::Edge;
  using Face_handle = typename VTr::Face_handle;
  using Context = typename VTr::Context;
  using Constraint_id = typename VTr::Constraint_id;

  // Validate polygons
  for (auto const& polygon : multipolygon.polygons_with_holes()) {
    if (!is_valid(polygon)) {
      std::cerr << "Invalid: polygon with holes" << std::endl;
      return false;
    }
  }

  std::vector<std::set<Constraint_id> > constraints;

  VTr vt;
  for (auto const& polygon : multipolygon.polygons_with_holes()) {
    std::set<Constraint_id> new_constraints;

    for (auto const& edge: polygon.outer_boundary().edges()) {
      try {
        Constraint_id cid = vt.insert_constraint(edge.source(), edge.target());
        new_constraints.insert(cid);
      } catch (typename VTr::Intersection_of_constraints_exception&) {
        std::cerr << "Invalid: (partly) overlapping polygons" << std::endl;
        return false;
      }
    }

    for (auto const& hole: polygon.holes()) {
      for (auto const& edge: hole.edges()) {
        try {
          Constraint_id cid = vt.insert_constraint(edge.source(), edge.target());
          new_constraints.insert(cid);
        } catch (typename VTr::Intersection_of_constraints_exception&) {
          std::cerr << "Invalid: hole boundary intersects something" << std::endl;
          return false;
        }
      }
    }

    constraints.emplace_back(std::move(new_constraints));
  }

  for (std::size_t i=0; i<constraints.size(); ++i) {
    for (auto const face: vt.all_face_handles()) {
      face->label() = 0;
      face->processed() = false;
    }

    struct Selected_constraints_checker
    {
      Selected_constraints_checker(const std::set<Constraint_id>& constraints)
        : m_constraints(constraints)
      { }

      bool operator()(const Edge& e, const VTr& vt) const
      {
        if (!vt.is_constrained(e)) {
          return false;
        }
        Vertex_handle va = e.first->vertex((e.second + 1)%3);
        Vertex_handle vb = e.first->vertex((e.second + 2)%3);
        for (const Context& c : vt.contexts(va, vb)) {
          if (m_constraints.count(c.id()) != 0) {
            return true;
          }
        }
        return false;
      };

      std::set<Constraint_id> m_constraints;
    };

    Selected_constraints_checker checker(constraints[i]);

    std::list<typename VTr::Face_handle> to_check;
    std::list<int> to_check_added_by;
    PR::label_region(vt, vt.infinite_face(), -1, to_check, to_check_added_by, checker); // exterior

    int regions = 0, holes = 0;
    while (!to_check.empty()) {
      if (to_check.front()->label() == 0) { // label = 0 means not labeled yet
        if (to_check_added_by.front() < 0) {
          PR::label_region(vt, to_check.front(), regions+1, to_check, to_check_added_by, checker);
          ++regions;
        } else {
          PR::label_region(vt, to_check.front(), -(holes+2), to_check, to_check_added_by, checker);
          ++holes;
        }
      }
      to_check.pop_front();
      to_check_added_by.pop_front();
    }

    // new constraints should be fully outside of the previous polygons
    for (const Edge& e : vt.all_edges()) {
      if (!vt.is_constrained(e)) {
        continue;
      }

      unsigned hits = 0;
      Vertex_handle va = e.first->vertex((e.second + 1)%3);
      Vertex_handle vb = e.first->vertex((e.second + 2)%3);
      for (const Context& c : vt.contexts(va, vb)) {
        if (constraints[i].count(c.id()) != 0) {
          ++hits;
        }
      }

      CGAL_assertion(hits == 0 || hits == 1);

      Face_handle fh = e.first;
      Face_handle nfh = e.first->neighbor(e.second);
      if (hits == 0) {
        if (fh->label() > 0 || nfh->label() > 0) {
          std::cerr << "Invalid: other polygon constraint touching selected constraint" << std::endl;
          return false;
        }
      } else if (hits == 1) {
        if (fh->label() > 0 && nfh->label() > 0) {
          std::cerr << "Invalid: polygon constraint not partly outside" << std::endl;
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
  using Validation_triangulation_base = CGAL::Constrained_triangulation_2<Kernel, Triangulation_data_structure, Validation_tag>;
  using Validation_triangulation = CGAL::Constrained_triangulation_plus_2<Validation_triangulation_base>;

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
  struct DefaultConstraintChecker {
    template <typename Edge, typename Tr>
    bool operator()(const Edge& e, const Tr& tt) const { return tt.is_constrained(e); }
  };

  template <class T,
            class ConstraintChecker = DefaultConstraintChecker>
  static void label_region(T& tt, Face_handle face, int label,
                           std::list<Face_handle>& to_check,
                           std::list<int>& to_check_added_by,
                           ConstraintChecker constraint_checker = ConstraintChecker{})
  {
    if (tt.dimension() < 2) {
      for (Face_handle fh : tt.all_face_handles()) {
        fh->label() = label;
      }
      return;
    }

    // std::cout << "Labeling region with " << label << std::endl;
    std::list<Face_handle> to_check_in_region;
    face->label() = label;
    to_check_in_region.push_back(face);
    face->processed() = true; // processed means added to a list (to ensure elements are only added once)

    while (!to_check_in_region.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!constraint_checker(Edge(to_check_in_region.front(), neighbour), tt)) {
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
      }
      to_check.pop_front();
      to_check_added_by.pop_front();
    }
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



/// \ingroup PkgPolygonRepairFunctions
/// computes the union of all polygons with holes in `p`
/// \tparam Kernel parameter of the input and output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
/// \pre Each polygon with holes must be free of self-intersections,
///      the outer boundaries must be counterclockwise and the holes clockwise oriented.
template <typename Kernel, typename Container>
Multipolygon_with_holes_2<Kernel,Container>
join(const Multipolygon_with_holes_2<Kernel,Container>& pa)
{
  struct Larger_than_zero {
    bool operator()(int i) const
    {
      return i > 0;
    }
  };

  CGAL::Polygon_repair::Boolean<Kernel, Container> bops;
  bops.insert(pa);
  bops.mark_domains();
  Larger_than_zero ltz;
  return bops(ltz);
}

/// \ingroup PkgPolygonRepairFunctions
/// computes the union of two polygonal domains
/// \tparam Kernel parameter of the output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
/// \tparam PA must be `Polygon_2<Kernel, Container>`, or `Polygon_with_holes_2<Kernel, Container>`, or `Multipolygon_with_holes_2<Kernel, Container>`
/// \tparam PB must be `Polygon_2<Kernel, Container>`, or `Polygon_with_holes_2<Kernel, Container>`, or `Multipolygon_with_holes_2<Kernel, Container>`
/// \pre The polygons `pa` and `pb` must be free of self-intersections,
///      the outer boundaries must be counterclockwise  and the holes clockwise oriented.
template <typename PA, typename PB, typename Kernel = Default, typename Container = Default>
#ifdef DOXYGEN_RUNNING
Multipolygon_with_holes_2<Kernel,Container>
#else
decltype(auto)
#endif
join(const PA& pa, const PB& pb, const Kernel& = Default(), const Container& = Default())
{
  typedef typename Default::Get<Kernel, typename PA::Traits>::type Traits;
  typedef typename Default::Get<Container, typename PA::Container>::type Container_;

  struct Larger_than_zero {
    bool operator()(int i) const
    {
      return i > 0;
    }
  };

  CGAL::Polygon_repair::Boolean<Traits,Container_> bops;
  bops.insert(pa);
  bops.insert(pb);
  bops.mark_domains();
  Larger_than_zero ltz;
  return bops(ltz);
}


/// \ingroup PkgPolygonRepairFunctions
/// computes the intersection of all polygons with holes in `p`
/// \tparam Kernel parameter of the input and output polygons. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
/// \pre Each polygon with holes must be free of self-intersections
///      the outer boundaries must be counterclockwise and the holes clockwise oriented.
template <typename Kernel, typename Container>
Multipolygon_with_holes_2<Kernel,Container>
intersect(const Multipolygon_with_holes_2<Kernel,Container>& p)
{
  struct Equal  {
    int val;
    Equal(int val)
    : val(val)
    {}

    bool operator()(int i) const
    {
      return i == val;
    }
  };

  CGAL::Polygon_repair::Boolean<Kernel,Container> bops;
  bops.insert(p);
  bops.mark_domains();
  Equal equal(p.number_of_polygons_with_holes());
  return bops(equal);
}


/// \ingroup PkgPolygonRepairFunctions
/// computes the intersection of two polygonal domains
/// \tparam Kernel parameter of the output polygon. Must be model of `ConstrainedDelaunayTriangulationTraits_2 `
/// \tparam Container parameter of the input and output polygons
/// \tparam PA must be `Polygon_2<Kernel, Container>`, or `Polygon_with_holes_2<Kernel, Container>`, or `Multipolygon_with_holes_2<Kernel, Container>`
/// \tparam PB must be `Polygon_2<Kernel, Container>`, or `Polygon_with_holes_2<Kernel, Container>`, or `Multipolygon_with_holes_2<Kernel, Container>`
/// \pre The polygons `pa` and `pb` must be free of self-intersections
///      the outer boundaries must be counterclockwise and the holes clockwise oriented.
template <typename PA, typename PB, typename Kernel = Default, typename Container = Default>
#ifdef DOXYGEN_RUNNING
Multipolygon_with_holes_2<Kernel, Container>
#else
decltype(auto)
#endif
intersect(const PA& pa, const PB& pb, const Kernel& = Default(), const Container& = Default())
{
  typedef typename Default::Get<Kernel, typename PA::Traits>::type Traits;
  typedef typename Default::Get<Container, typename PA::Container>::type Container_;

  struct Equal  {
    bool operator()(int i) const
    {
      return i == 2;
    }
  };

  CGAL::Polygon_repair::Boolean<Traits,Container_> bops;
  bops.insert(pa);
  bops.insert(pb);
  bops.mark_domains();
  Equal equal;
  return bops(equal);
}

} // namespace Polygon_repair
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_H
