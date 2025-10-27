// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
#define CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_triangulation_3/internal/config.h>

#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_data_3.h>
#include <CGAL/Triangulation_2/internal/Polyline_constraint_hierarchy_2.h>
#include <CGAL/Triangulation_segment_traverser_3.h>
#include <CGAL/unordered_flat_set.h>

#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/IO/File_binary_mesh_3.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/map.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <fstream>
#include <bitset>

#ifndef DOXYGEN_RUNNING

namespace CGAL {

template <typename T_3>
class Conforming_Delaunay_triangulation_3 : public T_3 {
public:
  static constexpr bool t_3_is_not_movable =
      CGAL::cdt_3_msvc_2019_or_older() || false == CGAL::is_nothrow_movable_v<T_3>;
  using Geom_traits = typename T_3::Geom_traits;
  using Vertex_handle = typename T_3::Vertex_handle;
  using Edge = typename T_3::Edge;
  using Facet = typename T_3::Facet;
  using Cell_handle = typename T_3::Cell_handle;
  using Point = typename T_3::Point;
  using Line = typename T_3::Geom_traits::Line_3;
  using Locate_type = typename T_3::Locate_type;

  inline static With_offset_tag with_offset{};
  inline static With_point_tag with_point{};
  inline static With_point_and_info_tag with_point_and_info{};

  Conforming_Delaunay_triangulation_3(const Geom_traits& gt = Geom_traits())
    : T_3(gt)
  {
  }

protected:
  struct Compare_vertex_handle {
    const T_3* tr;
    Compare_vertex_handle(const T_3* tr) : tr(tr) {}
    bool operator()(const Vertex_handle va, const Vertex_handle vb) const {
      return va < vb;
    }
  };

  using Constraint_hierarchy =
      Polyline_constraint_hierarchy_2<Vertex_handle, Compare_vertex_handle, Point>;
public:
  using Constrained_polyline_id = typename Constraint_hierarchy::Constraint_id;

protected:
  using Subconstraint = typename Constraint_hierarchy::Subconstraint;

  auto display_vert(Vertex_handle v) const{
    std::stringstream os;
    os.precision(17);
    os << IO::oformat(v, with_point);
    return os.str();
  }

  auto display_subcstr(Subconstraint subconstraint) const {
    auto [va, vb] = subconstraint;
    std::stringstream os;
    os << "(" << IO::oformat(va, with_offset) << ", " << IO::oformat(vb, with_offset) << ")"
       << ": [ " << display_vert(va) << " - " << display_vert(vb) << " ]";
    return os.str();
  }

  class Insert_in_conflict_visitor
  {
    Conforming_Delaunay_triangulation_3<T_3>* self;
  public:
    Insert_in_conflict_visitor(Conforming_Delaunay_triangulation_3* self) : self(self) {}

    template <class InputIterator>
    void process_cells_in_conflict(InputIterator cell_it, InputIterator end) {
      auto d = self->tr().dimension();
      for( ; cell_it != end; ++cell_it )
        for( int i = 0; i < d; ++i )
          for( int j = i+1; j <= d; ++j ) {
            auto v1 = (*cell_it)->vertex(i);
            auto v2 = (*cell_it)->vertex(j);
            if(self->tr().is_infinite(v1) || self->tr().is_infinite(v2)) continue;
            if(self->use_finite_edges_map()) {
              if(v1 > v2) std::swap(v1, v2);
              auto v1_index = v1->time_stamp();
              [[maybe_unused]] auto nb_erased = self->all_finite_edges[v1_index].erase(v2);
              if constexpr (cdt_3_can_use_cxx20_format()) if(self->debug_finite_edges_map() && nb_erased > 0) {
                std::cerr << cdt_3_format("erasing edge {} {}\n", self->display_vert((std::min)(v1, v2)),
                                        self->display_vert((std::max)(v1, v2)));
              }
            }
            auto [contexts_begin, contexts_end] =
              self->constraint_hierarchy.contexts(v1, v2);
            if(contexts_begin != contexts_end) {
              self->add_to_subconstraints_to_conform(v1, v2, contexts_begin->id());
            }
          }
    }

    void after_insertion(Vertex_handle v) const {
      if(!self->use_finite_edges_map()) return;
      CGAL_assertion(self->dimension() > 1);
      self->incident_edges(v, boost::make_function_output_iterator([this](Edge e) { self->new_edge(e); }));
      self->incident_cells(v, boost::make_function_output_iterator([this, v](Cell_handle c) {
        auto v_index = c->index(v);
        if(self->dimension() == 2) {
          auto j = self->cw(v_index);
          auto k = self->ccw(v_index);
          self->new_edge(Edge(c, j, k));
        } else {
          for(int i = 0; i < 3; ++i) {
            auto j = self->vertex_triple_index(v_index, i);
            auto k = self->vertex_triple_index(v_index, self->cw(i));
            self->new_edge(Edge(c, j, k));
          }
        }
      }));
    }

    void reinsert_vertices(Vertex_handle v) const { after_insertion(v); }
    Vertex_handle replace_vertex(Cell_handle c, int index, const Point& ) const
    {
      return c->vertex(index);
    }
    void hide_point(Cell_handle, const Point& ) const {}

    void insert_Steiner_point_on_constraint([[maybe_unused]] Constrained_polyline_id constraint,
                                            [[maybe_unused]] Vertex_handle va,
                                            [[maybe_unused]] Vertex_handle vb,
                                            [[maybe_unused]] Vertex_handle v_Steiner) const
    {
    }

    Vertex_handle insert_in_triangulation(const Point& p, Locate_type lt, Cell_handle c, int li, int lj) {
      return self->insert_impl_do_not_split(p, lt, c, li, lj, *this);
    }
  };

  auto make_subconstraint(Vertex_handle va, Vertex_handle vb) {
    return make_sorted_pair(va, vb);
  }

  void add_to_subconstraints_to_conform(Vertex_handle va, Vertex_handle vb,
                                        Constrained_polyline_id id) {
    const auto pair = make_subconstraint(va, vb);
#if CGAL_DEBUG_CDT_3 & 32
    std::cerr << "tr().subconstraints_to_conform.push("
              << display_subcstr(pair) << ")\n";
#endif // CGAL_DEBUG_CDT_3
    subconstraints_to_conform.push({pair, id});
  }

  template <typename Visitor>
  Constrained_polyline_id insert_constrained_edge_impl(Vertex_handle va, Vertex_handle vb,
                                             Visitor&) {
    if(va != vb) {
      if(segment_vertex_epsilon != 0.) {
        auto [min_dist, min_vertex] = min_distance_and_vertex_between_constraint_and_encroaching_vertex(va, vb);
        check_segment_vertex_distance_or_throw(va, vb, min_vertex, CGAL::to_double(min_dist),
                                               Check_distance::NON_SQUARED_DISTANCE);
      }
      const Constrained_polyline_id c_id = constraint_hierarchy.insert_constraint(va, vb);
      pair_of_vertices_to_cid.emplace(make_sorted_pair(va, vb), c_id);
      // traverse all the vertices along [va, vb] and add pairs of consecutive
      // vertices as sub-constraints.
      std::for_each(tr().segment_traverser_simplices_begin(va, vb), tr().segment_traverser_simplices_end(),
                    [&, prev = Vertex_handle{}](auto simplex) mutable {
                      // std::cerr << "- " << oformat(simplex, With_point_tag{}) << '\n';
                      if(simplex.dimension() == 0) {
                        const auto v = static_cast<Vertex_handle>(simplex);
                        v->ccdt_3_data().set_on_constraint(c_id);
                        if(prev != Vertex_handle{}) {
                          if(v != vb) {
                            v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_ON_EDGE);
                            constraint_hierarchy.add_Steiner(prev, vb, v);
                          }
                          add_to_subconstraints_to_conform(prev, v, c_id);
                        }
                        prev = v;
                      }
                    });
      return c_id;
    }
    return {};
  }

  void new_edge(Edge e)
  {
    if(!update_all_finite_edges_) return;
    auto [v1, v2] = make_sorted_pair(tr().vertices(e));
    if(tr().is_infinite(v1) || tr().is_infinite(v2))
      return;
    [[maybe_unused]] auto [_, inserted] = all_finite_edges[v1->time_stamp()].insert(v2);
    if constexpr (cdt_3_can_use_cxx20_format()) if (debug_finite_edges_map() && inserted) {
      if(v2 < v1) std::swap(v1, v2);
      std::cerr << cdt_3_format("new_edge({}, {})\n", display_vert(v1), display_vert(v2));
    }
  }

  void new_cell(Cell_handle c) {
    auto d = tr().dimension();
    for(int i = 0; i < d; ++i) {
      for(int j = i+1; j <= d; ++j) {
        new_edge(Edge(c, i, j));
      }
    }
  }

  auto new_cells_output_iterator()
  {
    return boost::function_output_iterator([this](Cell_handle c) {
      if(use_finite_edges_map())
        this->new_cell(c);
    });
  }

  template <typename Visitor>
  Vertex_handle insert_impl_do_not_split(const Point &p, Locate_type lt, Cell_handle c,
                                         int li, int lj, Visitor& visitor)
  {
    switch (tr().dimension()) {
    case 3: {
      typename T_3::Conflict_tester_3 tester(p, this);
      auto v = tr().insert_in_conflict(p, lt, c, li, lj, tester,
                                     visitor);
      new_vertex(v);
      return v;
    } // dim 3
    case 2: {
      typename T_3::Conflict_tester_2 tester(p, this);
      auto v = tr().insert_in_conflict(p, lt, c, li, lj, tester, visitor);
      if(use_finite_edges_map()) {
        new_vertex(v);
        tr().incident_edges(v, boost::make_function_output_iterator([&](Edge e) { this->new_edge(e); }));
      }
      return v;
    } // dim 2
    default:
      // dimension <= 1
      // Do not use the generic insert.
      auto v = tr().insert(p, c);
      if(use_finite_edges_map()) {
        new_vertex(v);
        all_finite_edges.clear();
        if (debug_finite_edges_map()) std::cerr << "all_finite_edges.clear()\n";
        for(auto e: tr().all_edges()) {
          new_edge(e);
        }
      }
      return v;
    }
  }

  template <typename Visitor>
  Vertex_handle insert_impl(const Point &p, Locate_type lt, Cell_handle c,
                            int li, int lj, Visitor& visitor)
  {
    Vertex_handle v1, v2;
    bool split_constrained_edge = false;
    if(lt == T_3::EDGE) {
      v1 = c->vertex(li);
      v2 = c->vertex(lj);
      if(constraint_hierarchy.is_subconstraint(v1, v2)) {
        split_constrained_edge = true;
      }
    }
    Vertex_handle new_vertex = insert_impl_do_not_split(p, lt, c, li, lj, visitor);
    if(split_constrained_edge) {
      constraint_hierarchy.split_constraint(v1, v2, new_vertex);
    }
    return new_vertex;
  }

  void update_bbox(const Point& p) {
    bbox = bbox + p.bbox();
    if(max_bbox_edge_length) {
      update_max_bbox_edge_length();
    }
  }

  void update_max_bbox_edge_length() {
    double d_x = bbox.xmax() - bbox.xmin();
    double d_y = bbox.ymax() - bbox.ymin();
    double d_z = bbox.zmax() - bbox.zmin();

    max_bbox_edge_length = (std::max)(d_x, (std::max)(d_y, d_z));
  }

public:
  void set_segment_vertex_epsilon(double epsilon) {
    segment_vertex_epsilon = epsilon;
  }

  bool debug_Steiner_points() const {
    return debug_flags[static_cast<int>(Debug_flags::Steiner_points)];
  }

  void debug_Steiner_points(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::Steiner_points), b);
  }

  bool debug_input_faces() const {
    return debug_flags[static_cast<int>(Debug_flags::input_faces)];
  }

  void debug_input_faces(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::input_faces), b);
  }

  bool debug_missing_region() const {
    return debug_flags[static_cast<int>(Debug_flags::missing_region)];
  }

  void debug_missing_region(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::missing_region), b);
  }

  bool debug_regions() const {
    return debug_flags[static_cast<int>(Debug_flags::regions)];
  }

  void debug_regions(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::regions), b);
  }

  bool debug_copy_triangulation_into_hole() const {
    return debug_flags[static_cast<int>(Debug_flags::copy_triangulation_into_hole)];
  }

  void debug_copy_triangulation_into_hole(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::copy_triangulation_into_hole), b);
  }

  bool debug_validity() const {
    return debug_flags[static_cast<int>(Debug_flags::validity)];
  }

  void debug_validity(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::validity), b);
  }

  bool use_older_cavity_algorithm() const {
    return debug_flags[static_cast<int>(Debug_flags::use_older_cavity_algorithm)];
  }

  void use_older_cavity_algorithm(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::use_older_cavity_algorithm), b);
  }

  bool debug_finite_edges_map() const {
    return debug_flags[static_cast<int>(Debug_flags::debug_finite_edges_map)];
  }

  void debug_finite_edges_map(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::debug_finite_edges_map), b);
  }

  bool use_finite_edges_map() const {
    return update_all_finite_edges_ && debug_flags[static_cast<int>(Debug_flags::use_finite_edges_map)];
  }

  void use_finite_edges_map(bool b) {
    debug_flags.set(static_cast<int>(Debug_flags::use_finite_edges_map), b);
  }

  Vertex_handle insert(const Point &p, Locate_type lt, Cell_handle c,
                       int li, int lj)
  {
    update_bbox(p);
    auto v = insert_impl(p, lt, c, li, lj, insert_in_conflict_visitor);
    restore_Delaunay(insert_in_conflict_visitor);
    return v;
  }

  Vertex_handle insert(const Point &p, Cell_handle start = {}) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr().locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }

  Constrained_polyline_id insert_constrained_edge(Vertex_handle va, Vertex_handle vb)
  {
    const auto id = insert_constrained_edge_impl(va, vb, insert_in_conflict_visitor);
    restore_Delaunay(insert_in_conflict_visitor);
    return id;
  }

  bool is_edge(Vertex_handle va, Vertex_handle vb) const {
    const bool is_edge_v1 =
        ((debug_finite_edges_map() && use_finite_edges_map()) || !use_finite_edges_map()) && tr().tds().is_edge(va, vb);

    if(use_finite_edges_map() && va > vb) std::swap(va, vb);
    const auto va_index = va->time_stamp();
    const bool is_edge_v2 =
        use_finite_edges_map() && all_finite_edges[va_index].find(vb) != all_finite_edges[va_index].end();

    if(debug_finite_edges_map() && use_finite_edges_map() && is_edge_v1 != is_edge_v2) {
      std::cerr << "!! Inconsistent edge status\n";
      std::cerr << "  -> constraint " << display_vert(va) << "     " << display_vert(vb) << '\n';
      std::cerr << "  ->     edge " << (is_edge_v1 ? "is" : "is not") << " in the triangulation\n";
      std::cerr << "  ->     edge " << (is_edge_v2 ? "is" : "is not") << " in all_finite_edges\n";
      debug_dump("bug-inconsistent-edge-status");
      CGAL_error();
    }
    const bool is_edge = use_finite_edges_map() ? is_edge_v2 : is_edge_v1;
    return is_edge;
  }

  using T_3::is_edge;

  bool is_conforming() const {
    return std::all_of(constraint_hierarchy.subconstraints_begin(),
                       constraint_hierarchy.subconstraints_end(),
                       [this](const auto &sc) {
                         const auto [va, vb] = sc;
                         const auto is_edge = this->is_edge(va, vb);
#if CGAL_DEBUG_CDT_3 & 128 && CGAL_CAN_USE_CXX20_FORMAT
                         std::cerr << cdt_3_format("is_conforming>> Edge is 3D: {}  ({} , {})\n",
                                                  is_edge,
                                                  CGAL::IO::oformat(va, with_point_and_info),
                                                  CGAL::IO::oformat(vb, with_point_and_info));
#endif // CGAL_DEBUG_CDT_3
                         return is_edge;
                       });
  }

  enum class Check_distance { SQUARED_DISTANCE, NON_SQUARED_DISTANCE };

  void check_segment_vertex_distance_or_throw(Vertex_handle va,
                                              Vertex_handle vb,
                                              Vertex_handle min_vertex,
                                              double min_dist,
                                              Check_distance option)
  {
    if(!max_bbox_edge_length) {
      update_max_bbox_edge_length();
    }
    if((option == Check_distance::NON_SQUARED_DISTANCE && min_dist < segment_vertex_epsilon * *max_bbox_edge_length) ||
       (option == Check_distance::SQUARED_DISTANCE &&
        min_dist < CGAL::square(segment_vertex_epsilon * *max_bbox_edge_length)))
    {
      std::stringstream ss;
      ss.precision(std::cerr.precision());
      ss << "A constrained segment is too close to a vertex.\n";
      ss << "  -> vertex " << display_vert(min_vertex) << '\n';
      ss << "  -> constrained segment " << display_vert(va) << "  -  " << display_vert(vb) << '\n';
      ss << "  -> distance = " << min_dist << '\n';
      ss << "  -> max_bbox_edge_length = " << *max_bbox_edge_length << '\n';
      CGAL_error_msg(ss.str().c_str());
    }
  }

  auto ancestors_of_Steiner_vertex_on_edge(Vertex_handle v) const {
    std::pair<Vertex_handle, Vertex_handle> result;
    CGAL_precondition(v->ccdt_3_data().is_Steiner_vertex_on_edge());
    CGAL_assertion(v->ccdt_3_data().number_of_incident_constraints() == 1);
    const auto v_time_stamp = v->time_stamp();
    const auto constrained_polyline_id = v->ccdt_3_data().constrained_polyline_id(*this);
    const auto first = this->constraint_hierarchy.vertices_in_constraint_begin(constrained_polyline_id);
    const auto end =   this->constraint_hierarchy.  vertices_in_constraint_end(constrained_polyline_id);
    std::cerr << "ancestors_of_Steiner_vertex_on_edge " << display_vert(v) << '\n';
    for(auto it = first; it != end; ++it) {
      std::cerr << "  - " << display_vert(*it) << '\n';
    }
    CGAL_assertion(first != end);
    const auto last = std::prev(this->constraint_hierarchy.vertices_in_constraint_end(constrained_polyline_id));
    CGAL_assertion(first != last);
    CGAL_assertion((*first)->time_stamp() < v_time_stamp);
    CGAL_assertion((*last)->time_stamp() < v_time_stamp);
    for(auto it = first; (*it) != v; ++it) {
      if((*it)->time_stamp() < v_time_stamp) {
        result.first = *it;
      }
      CGAL_assertion(it != last);
    }
    for(auto it = last; (*it) != v; --it) {
      if((*it)->time_stamp() < v_time_stamp) {
        result.second = *it;
      }
      CGAL_assertion(it != first);
    }
    return result;
  }

  auto min_distance_and_vertex_between_constraint_and_encroaching_vertex(Vertex_handle va, Vertex_handle vb) const {
    struct Result {
      typename T_3::Geom_traits::FT min_dist = (std::numeric_limits<double>::max)();
      Vertex_handle v = {};
    } result;
    const auto vector_of_encroaching_vertices = encroaching_vertices(va, vb);
    for(auto v: vector_of_encroaching_vertices) {
      const auto dist = CGAL::approximate_sqrt(squared_distance(tr().point(v), Line{tr().point(va), tr().point(vb)}));
      if(dist < result.min_dist) {
        result = Result{dist, v};
      }
    }
    return result;
  }

  bool write_missing_segments_file(std::ostream &out) {
    bool any_missing_segment = false;
    std::for_each(
        constraint_hierarchy.subconstraints_begin(), constraint_hierarchy.subconstraints_end(),
        [this, &out, &any_missing_segment](const auto &sc) {
          const auto [v0, v1] = sc;
          if (!this->is_edge(v0, v1)) {
            out << "2 " << this->tr().point(v0) << " " << this->tr().point(v1)
                << '\n';
            any_missing_segment = true;
          }
        });
    return any_missing_segment;
  }

  void write_all_segments_file(std::ostream &out) {
    std::for_each(
        constraint_hierarchy.subconstraints_begin(), constraint_hierarchy.subconstraints_end(),
        [this, &out](const auto &sc) {
          const auto [v0, v1] = sc;
          out << "2 " << this->tr().point(v0) << " " << this->tr().point(v1) << '\n';
        });
  }

  /// @{
  /// remove functions cannot be called
  void remove(Vertex_handle) = delete;
  void remove_cluster() = delete;
  /// @}

  static std::string io_signature() {
    return Get_io_signature<T_3>()();
  }
protected:
  void debug_dump(std::string filename) const {
    {
      std::ofstream dump(filename + ".binary.cgal", std::ios::binary);
      CGAL::IO::save_binary_file(dump, *this);
    }
    {
      std::ofstream dump(filename + "_point.xyz");
      dump.precision(17);
      for(auto vh: this->finite_vertex_handles()){
        dump << this->point(vh) << '\n';
      }
    }
  }

  template <typename Visitor>
  void restore_Delaunay(Visitor& visitor) {
    update_all_finite_edges();
    while(!subconstraints_to_conform.empty()) {
      const auto [subconstraint, constrained_polyline_id] = subconstraints_to_conform.top();
      subconstraints_to_conform.pop();
      const auto [va, vb] = subconstraint;
      if(!constraint_hierarchy.is_subconstraint(va, vb)) {
        continue;
      }
#if CGAL_DEBUG_CDT_3 & 32
      std::cerr << "tr().subconstraints_to_conform.pop()="
                << display_subcstr(subconstraint) << "\n";
#endif // CGAL_DEBUG_CDT_3
      conform_subconstraint(subconstraint, constrained_polyline_id, visitor);
    }
  }

  template <typename Visitor>
  Vertex_handle insert_Steiner_point_on_subconstraint(
      Point steiner_pt, Cell_handle hint,
      Subconstraint subconstraint, Constrained_polyline_id constraint, Visitor& visitor)
  {
    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    Locate_type lt;
    int li, lj;
    const Cell_handle c = tr().locate(steiner_pt, lt, li, lj, hint);
    const Vertex_handle v = visitor.insert_in_triangulation(steiner_pt, lt, c, li, lj);
    v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_ON_EDGE);
    if(lt != T_3::VERTEX) {
      v->ccdt_3_data().set_on_constraint(constraint);
    }
    constraint_hierarchy.add_Steiner(va, vb, v);
    visitor.insert_Steiner_point_on_constraint(constraint, va, vb, v);
    add_to_subconstraints_to_conform(va, v, constraint);
    add_to_subconstraints_to_conform(v, vb, constraint);

    new_vertex(v);

    return v;
  }

  /// Return `true` if a Steiner point was inserted
  template <typename Visitor>
  bool conform_subconstraint(Subconstraint subconstraint,
                             Constrained_polyline_id constraint,
                             Visitor& visitor)
  {
    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    CGAL_assertion(va != vb);
    if(!this->is_edge(va, vb)) {
      const auto& [steiner_pt, hint, ref_vertex] = construct_Steiner_point(constraint, subconstraint);
      [[maybe_unused]] const auto v =
          insert_Steiner_point_on_subconstraint(steiner_pt, hint, subconstraint, constraint, visitor);
      if(debug_Steiner_points()) {
        const auto [c_start, c_end] = constraint_extremities(constraint);
        std::cerr << "(" << IO::oformat(va, with_offset) << ", " << IO::oformat(vb, with_offset) << ")";
        std::cerr << ": [ " << display_vert(c_start) << " - " << display_vert(c_end) << " ] ";
        std::cerr << "  new vertex " << display_vert(v) << '\n';
      }
      return true;
    }

    return false;
  }

  Constrained_polyline_id constraint_from_extremities(Vertex_handle va, Vertex_handle vb) const {
    if(va->ccdt_3_data().number_of_incident_constraints() == 0 || vb->ccdt_3_data().number_of_incident_constraints() == 0)
    {
      return {};
    }
    auto it = pair_of_vertices_to_cid.find(make_sorted_pair(va, vb));
    if(it != pair_of_vertices_to_cid.end()) {
      return it->second;
    }
    return {};
    // @TODO: cleanup the rest of the function, and `constraint_around`
    Constrained_polyline_id c_id = constraint_around(va, vb, false);
    if(c_id != Constrained_polyline_id{}) return c_id;
    c_id = constraint_around(vb, va, false);
    if(c_id != Constrained_polyline_id{}) return c_id;
    c_id = constraint_around(va, vb, true);
    return c_id;
  }

  auto constraint_extremities(Constrained_polyline_id c_id) const {
      CGAL_assertion(std::find(this->constraint_hierarchy.constraints_begin(),
                               this->constraint_hierarchy.constraints_end(), c_id) != this->constraint_hierarchy.constraints_end());
      CGAL_assertion(this->constraint_hierarchy.vertices_in_constraint_begin(c_id) !=
                     this->constraint_hierarchy.vertices_in_constraint_end(c_id));
#if CGAL_DEBUG_CDT_3 & 8
      std::cerr << "constraint " << (void*) c_id.vl_ptr() << " has "
                << c_id.vl_ptr()->skip_size() << " vertices\n";
#endif // CGAL_DEBUG_CDT_3
      const auto begin = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
      const auto end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
      const auto c_va = *begin;
      const auto c_vb = *std::prev(end);
    return std::make_pair(c_va, c_vb);
  }

  Constrained_polyline_id constraint_around(Vertex_handle va, Vertex_handle vb, bool expensive = true) const {
    auto constraint_id_goes_to_vb = [this, va, vb](Constrained_polyline_id c_id) {
      const auto [c_va, c_vb] = constraint_extremities(c_id);
      if (va == c_va && vb == c_vb)
        return true;
      if (vb == c_va && va == c_vb)
        return true;
      return false;
    }; // end lambda constraint_id_goes_to_vb
    if (va->ccdt_3_data().number_of_incident_constraints() == 1)
    {
      const Constrained_polyline_id c_id = va->ccdt_3_data().constrained_polyline_id(*this);
      CGAL_assertion(c_id != Constrained_polyline_id{});
      if(constraint_id_goes_to_vb(c_id)) return c_id;
    } else if (expensive == true && va->ccdt_3_data().number_of_incident_constraints() > 1) {
      boost::container::small_vector<Vertex_handle, 64> adj_vertices;
      this->finite_adjacent_vertices(va, std::back_inserter(adj_vertices));
      for(auto other_v: adj_vertices) {
        for(auto context: this->constraint_hierarchy.contexts(va, other_v)) {
          const Constrained_polyline_id c_id = context.id();
          if(constraint_id_goes_to_vb(c_id)) return c_id;
        }
      }
    }
    return Constrained_polyline_id{};
  }

  struct Construct_Steiner_point_return_type {
    typename T_3::Geom_traits::Point_3 point;
    Cell_handle hint;
    Vertex_handle reference_vertex;
  };

  auto encroaching_vertices(Vertex_handle va, Vertex_handle vb) const {
    auto& gt = tr().geom_traits();
    auto angle_functor = gt.angle_3_object();

    const auto& pa = tr().point(va);
    const auto& pb = tr().point(vb);

    namespace bc = boost::container;
    bc::flat_set<Vertex_handle, std::less<Vertex_handle>,
                 bc::small_vector<Vertex_handle, 256>>
        encroaching_vertices;
    auto register_vertex = [this,&encroaching_vertices](Vertex_handle v) {
      if(tr().is_infinite(v)) return;
      // std::cerr << "register_vertex " << display_vert(v) << '\n';
      encroaching_vertices.insert(v);
    };
    auto fill_encroaching_vertices = [&](const auto simplex) {
#if CGAL_DEBUG_CDT_3 & 0x10
      std::cerr << " - " << IO::oformat(simplex, With_point_tag{}) << '\n';
#endif // CGAL_DEBUG_CDT_3
      auto visit_cell = [&](Cell_handle cell) {
        for(int i = 0, end = this->tr().dimension() + 1; i < end; ++i) {
          const auto v = cell->vertex(i);
          register_vertex(v);
        }
      };
      switch(simplex.dimension()) {
      case 3: {
        const auto cell = static_cast<Cell_handle>(simplex);
        visit_cell(cell);
      } break;
      case 2: {
        const auto [cell, facet_index] = static_cast<Facet>(simplex);
        visit_cell(cell);
        if(tr().dimension() > 2) {
          const auto [other_cell, other_index] = tr().mirror_facet({cell, facet_index});
          register_vertex(other_cell->vertex(other_index));
        }
        break;
      }
      case 1: {
        auto edge = static_cast<Edge>(simplex);
        if(tr().dimension() < 3) {
          auto [cell, i, j] = edge;
          visit_cell(cell);
          if(tr().dimension() < 2) break;
          auto neighbor_cell = cell->neighbor(3 - i - j);
          visit_cell(neighbor_cell);
          break;
        }
        auto circ = tr().incident_cells(edge);
        CGAL_assertion(circ != nullptr);
        const auto end = circ;
        do {
          visit_cell(circ);
        } while(++circ != end);
      } break;
      case 0: {
        const auto v = static_cast<Vertex_handle>(simplex);
        if(v != va && v != vb) {
          std::cerr << "!! The constraint passes through a vertex!\n";
          std::cerr << "  -> constraint " << display_vert(va) << "     " << display_vert(vb) << '\n';
          std::cerr << "  ->     vertex " << display_vert(v) << '\n';
#if CGAL_DEBUG_CDT_3
          debug_dump("bug-through-vertex");
#endif
          CGAL_error();
        }
      } break;
      default: CGAL_unreachable();
      } // end switch
    };
    std::for_each(tr().segment_traverser_simplices_begin(va, vb), tr().segment_traverser_simplices_end(),
                  fill_encroaching_vertices);
    auto vector_of_encroaching_vertices = encroaching_vertices.extract_sequence();
#if CGAL_DEBUG_CDT_3 & 0x10
    std::cerr << "  -> vector_of_encroaching_vertices (before filter):\n";
    std::for_each(vector_of_encroaching_vertices.begin(),
                  vector_of_encroaching_vertices.end(),
                  [this](Vertex_handle v){
                    std::cerr << "    " << this->display_vert(v) << '\n';
                  });
#endif // CGAL_DEBUG_CDT_3
    auto end = std::remove_if(vector_of_encroaching_vertices.begin(),
                              vector_of_encroaching_vertices.end(),
                              [va, vb, pa, pb, &angle_functor, this](Vertex_handle v) {
                                if(va == v || vb == v) return true;
                                return angle_functor(pa,
                                                    this->tr().point(v),
                                                    pb) == ACUTE;
                              });
#if CGAL_DEBUG_CDT_3 & 0x10
    std::cerr << "  -> vector_of_encroaching_vertices (after filter):\n";
    std::for_each(vector_of_encroaching_vertices.begin(), end, [&](Vertex_handle v) {
      std::cerr << "    " << this->display_vert(v) << "  angle " << approximate_angle(pa, this->tr().point(v), pb)
                << '\n';
    });
#endif // CGAL_DEBUG_CDT_3
    vector_of_encroaching_vertices.erase(end, vector_of_encroaching_vertices.end());
    return vector_of_encroaching_vertices;
  }

  Construct_Steiner_point_return_type
  construct_Steiner_point(Constrained_polyline_id constrained_polyline_id, Subconstraint subconstraint)
  {
    auto& gt = tr().geom_traits();
    auto compare_angle_functor = gt.compare_angle_3_object();
    auto vector_functor = gt.construct_vector_3_object();
    auto midpoint_functor = gt.construct_midpoint_3_object();
    auto scaled_vector_functor = gt.construct_scaled_vector_3_object();
    auto sq_length_functor = gt.compute_squared_length_3_object();
    auto sc_product_functor = gt.compute_scalar_product_3_object();
    auto translate_functor = gt.construct_translated_point_3_object();

    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    const auto& pa = tr().point(va);
    const auto& pb = tr().point(vb);
    const auto [orig_va, orig_vb] = constraint_extremities(constrained_polyline_id);
    const auto& orig_pa = tr().point(orig_va);
    const auto& orig_pb = tr().point(orig_vb);

    if(this->dimension() < 2) {
      std::cerr << "dim < 2: midpoint\n";
      return {midpoint_functor(pa, pb), va->cell(), va};
    }

#if CGAL_DEBUG_CDT_3 & 0x10
    std::cerr << "construct_Steiner_point( " << display_vert(va) << " , "
              << display_vert(vb) << " )\n";
#endif // CGAL_DEBUG_CDT_3

    const auto vector_of_encroaching_vertices = encroaching_vertices(va, vb);
    CGAL_assertion(vector_of_encroaching_vertices.size() > 0);

    const auto reference_vertex_it = std::max_element(
        vector_of_encroaching_vertices.begin(), vector_of_encroaching_vertices.end(),
        [pa, pb, &compare_angle_functor, this](Vertex_handle v1,
                                               Vertex_handle v2) {
          return compare_angle_functor(pa, this->tr().point(v1), pb,
                                       pa, this->tr().point(v2), pb) == SMALLER;
        });
    CGAL_assertion(reference_vertex_it != vector_of_encroaching_vertices.end());
#if CGAL_CDT_3_DEBUG_CONFORMING
    std::cerr << "  -> reference point: " << display_vert(*reference_vertex_it)
              << '\n';
#endif // CGAL_CDT_3_DEBUG_CONFORMING
    const auto reference_vertex = *reference_vertex_it;
    const auto& reference_point = tr().point(reference_vertex);

    const auto vector_ab = vector_functor(pa, pb);

    if(reference_vertex->ccdt_3_data().is_Steiner_vertex_on_edge()) {
      CGAL_assertion(reference_vertex->ccdt_3_data().number_of_incident_constraints() == 1);
      const auto ref_constrained_polyline_id = reference_vertex->ccdt_3_data().constrained_polyline_id(*this);
      const auto [ref_va, ref_vb] = constraint_extremities(ref_constrained_polyline_id);
#if CGAL_CDT_3_DEBUG_CONFORMING
      std::cerr << "  reference point is on constraint: " << display_vert(ref_va)
                << "    " << display_vert(ref_vb) << '\n'
                << "  original constraint:              " << display_vert(orig_va)
                << "    " << display_vert(orig_vb) << '\n';
#endif // CGAL_CDT_3_DEBUG_CONFORMING
      const auto vector_orig_ab = vector_functor(orig_pa, orig_pb);
      const auto length_ab = CGAL::approximate_sqrt(sq_length_functor(vector_ab));
      auto return_orig_result_point =
          [&](auto lambda, Point orig_pa, Point orig_pb)
              -> Construct_Steiner_point_return_type
          {
            const auto vector_orig_ab = vector_functor(orig_pa, orig_pb);
            const auto inter_point = translate_functor(orig_pa, scaled_vector_functor(vector_orig_ab, lambda));
            const auto dist_a_result = CGAL::approximate_sqrt(sq_length_functor(vector_functor(pa, inter_point)));
            const auto ratio = dist_a_result / length_ab;
            const auto result_point = (ratio < 0.2 || ratio > 0.8)
                                          ? midpoint_functor(pa, pb)
                                          : inter_point;

#if CGAL_CDT_3_DEBUG_CONFORMING
            std::cerr << "  ref ratio = " << ratio << '\n';
            std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_CDT_3_DEBUG_CONFORMING
            return {result_point, reference_vertex->cell(), reference_vertex};
          };

      const auto length_orig_ab = CGAL::approximate_sqrt(sq_length_functor(vector_orig_ab));
      if(ref_va == orig_va || ref_vb == orig_va) {
        const auto vector_orig_a_ref = vector_functor(orig_pa, reference_point);
        const auto length_orig_a_ref = CGAL::approximate_sqrt(sq_length_functor(vector_orig_a_ref));
        const auto lambda = length_orig_a_ref / length_orig_ab;
        return return_orig_result_point(lambda, orig_pa, orig_pb);
      } else if(ref_va == orig_vb || ref_vb == orig_vb) {
        const auto vector_orig_b_ref = vector_functor(orig_pb, reference_point);
        const auto length_orig_b_ref = CGAL::approximate_sqrt(sq_length_functor(vector_orig_b_ref));
        const auto lambda = length_orig_b_ref / length_orig_ab;
        return return_orig_result_point(lambda, orig_pb, orig_pa);
      }
    } else {
      if(segment_vertex_epsilon > 0) {
        if(!max_bbox_edge_length) {
          update_max_bbox_edge_length();
        }
        auto sq_dist = squared_distance(reference_point, Line{orig_pa, orig_pb});
        check_segment_vertex_distance_or_throw(orig_va, orig_vb, reference_vertex, CGAL::to_double(sq_dist),
                                               Check_distance::SQUARED_DISTANCE);
      }
    }
    // compute the projection of the reference point
    const auto vector_a_ref = vector_functor(pa, reference_point);
    const auto lambda = sc_product_functor(vector_a_ref, vector_ab) / sq_length_functor(vector_ab);
    const auto result_point = (lambda < 0.2 || lambda > 0.8)
                                  ? midpoint_functor(pa, pb)
                                  : translate_functor(pa, scaled_vector_functor(vector_ab, lambda));

#if CGAL_CDT_3_DEBUG_CONFORMING
    std::cerr << "  lambda = " << lambda << '\n';
    std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_CDT_3_DEBUG_CONFORMING
    return {result_point, reference_vertex->cell(), reference_vertex};
  }

protected:
  T_3& tr() { return *this; };
  const T_3& tr() const { return *this; };

  Compare_vertex_handle comp = {this};
  Constraint_hierarchy constraint_hierarchy = {comp};
  static_assert(CGAL::cdt_3_msvc_2019_or_older() || CGAL::is_nothrow_movable_v<Constraint_hierarchy>);
  Bbox_3 bbox{};
  double segment_vertex_epsilon = 1e-8;
  std::optional<double> max_bbox_edge_length;
  using Pair_of_vertex_handles = std::pair<Vertex_handle, Vertex_handle>;
  boost::container::map<Pair_of_vertex_handles, Constrained_polyline_id> pair_of_vertices_to_cid;
  Insert_in_conflict_visitor insert_in_conflict_visitor = {this};

  using Stack_info = std::pair<Subconstraint, Constrained_polyline_id>;
  using Subconstraints_to_conform = std::stack<Stack_info, std::vector<Stack_info>>;
  Subconstraints_to_conform subconstraints_to_conform;

  std::vector<CGAL::unordered_flat_set<Vertex_handle>> all_finite_edges;
  bool update_all_finite_edges_ = false;

  void update_all_finite_edges() {
    if(!update_all_finite_edges_) {
      update_all_finite_edges_ = true;
      if(use_finite_edges_map()) {
        all_finite_edges.clear();
        all_finite_edges.resize(tr().number_of_vertices()+1);
        for(auto e: tr().all_edges()) {
          new_edge(e);
        }
      }
    }
  }

  void new_vertex(Vertex_handle v) {
    if(use_finite_edges_map() && v->time_stamp() >= all_finite_edges.size()) {
      all_finite_edges.emplace_back();
      CGAL_assertion(v->time_stamp() == all_finite_edges.size() - 1);
    }
  }

  enum class Debug_flags {
    Steiner_points = 0,
    conforming,
    input_faces,
    missing_region,
    regions,
    copy_triangulation_into_hole,
    validity,
    use_older_cavity_algorithm,
    debug_finite_edges_map,
    use_finite_edges_map,
    nb_of_flags
  };
  std::bitset<static_cast<int>(Debug_flags::nb_of_flags)> debug_flags{};
  bool is_Delaunay = true;
};

} // end CGAL

#endif // not DOXYGEN_RUNNING

#

#endif // CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
