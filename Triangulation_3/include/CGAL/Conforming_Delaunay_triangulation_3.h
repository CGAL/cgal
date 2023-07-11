// Copyright (c) 2019  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
#define CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/Triangulation_2/internal/Polyline_constraint_hierarchy_2.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/IO/File_binary_mesh_3.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp> /// @TODO Requires Boost 1.66

#include <fstream>

namespace CGAL {

enum class CDT_3_vertex_type { FREE, CORNER, STEINER_ON_EDGE, STEINER_IN_FACE };

template <typename Gt, typename Vb = Triangulation_vertex_base_3<Gt> >
class Conforming_Delaunay_triangulation_vertex_base_3 : public Vb {
public:
  int nb_of_incident_constraints = 0;
  void* c_id = nullptr;

  // To get correct vertex type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Conforming_Delaunay_triangulation_vertex_base_3 <Gt, Vb3> Other;
  };

  using Vb::Vb;

  template<typename Triangulation>
  auto constraint_id(const Triangulation&) const {
    using C_id = typename Triangulation::Constraint_id;
    using Vertex_list_ptr = decltype(std::declval<C_id>().vl_ptr());
    return C_id(static_cast<Vertex_list_ptr>(c_id));
  }

  CDT_3_vertex_type vertex_type() const { return m_vertex_type; }
  void set_vertex_type(CDT_3_vertex_type type) { m_vertex_type = type; }
  bool is_Steiner_vertex_on_edge() const { return m_vertex_type == CDT_3_vertex_type::STEINER_ON_EDGE; }
  bool is_Steiner_vertex_in_face() const { return m_vertex_type == CDT_3_vertex_type::STEINER_IN_FACE; }

  static std::string io_signature() {
    return Get_io_signature<Vb>()();
  }
private:
  CDT_3_vertex_type m_vertex_type = CDT_3_vertex_type::FREE;
};


template <typename T_3>
class Conforming_Delaunay_triangulation_3 : public T_3 {
public:
  using Vertex_handle = typename T_3::Vertex_handle;
  using Edge = typename T_3::Edge;
  using Facet = typename T_3::Facet;
  using Cell_handle = typename T_3::Cell_handle;
  using Point = typename T_3::Point;
  using Locate_type = typename T_3::Locate_type;

protected:
  struct Compare_vertex_handle {
    const T_3* tr;
    Compare_vertex_handle(const T_3* tr) : tr(tr) {}
    bool operator()(const Vertex_handle va, const Vertex_handle vb) const {
      return tr->compare_xyz(tr->point(va), tr->point(vb)) == SMALLER;
    }
  };

  using Constraint_hierarchy =
      Polyline_constraint_hierarchy_2<Vertex_handle, Compare_vertex_handle, const Point&>;
public:
  using Constraint_id = typename Constraint_hierarchy::Constraint_id;

protected:
  using Subconstraint = typename Constraint_hierarchy::Subconstraint;

#if CGAL_DEBUG_CDT_3
  auto display_vert(Vertex_handle v) {
    std::stringstream os;
    os.precision(17);
    os << IO::oformat(v) << "=(" << this->tr.point(v) << ")";
    return os.str();
  };
  auto display_subcstr(Subconstraint subconstraint) {
    std::stringstream os;
    os << "[ " << display_vert(subconstraint.first)
       << " - " << display_vert(subconstraint.second) << " ]";
    return os.str();
  };
#endif // CGAL_DEBUG_CDT_3

  class Insert_in_conflict_visitor
  {
    Conforming_Delaunay_triangulation_3<T_3>& self;
  public:
    Insert_in_conflict_visitor(Conforming_Delaunay_triangulation_3& self) : self(self) {}

    template <class InputIterator>
    void process_cells_in_conflict(InputIterator cell_it, InputIterator end) {
      for( ; cell_it != end; ++cell_it )
        for( int i = 0; i < self.tr.dimension(); ++i )
          for( int j = i+1; j < self.tr.dimension()+1; ++j ) {
            auto v1 = (*cell_it)->vertex(i);
            auto v2 = (*cell_it)->vertex(j);
            if(self.tr.is_infinite(v1) || self.tr.is_infinite(v2)) continue;
            auto [contexts_begin, contexts_end] =
              self.constraint_hierarchy.contexts_range(v1, v2);
            if(contexts_begin != contexts_end) {
              self.add_to_subconstraints_to_conform(v1, v2,
                                                    contexts_begin->id());
            }
          }
    }
    void reinsert_vertices(Vertex_handle) const {}
    Vertex_handle replace_vertex(Cell_handle c, int index, const Point& ) const
    {
      return c->vertex(index);
    }
    void hide_point(Cell_handle, const Point& ) const {}

    void insert_Steiner_point_on_constraint([[maybe_unused]] Constraint_id constraint,
                                            [[maybe_unused]] Vertex_handle va,
                                            [[maybe_unused]] Vertex_handle vb,
                                            [[maybe_unused]] Vertex_handle v_Steiner) const
    {
    }
  };

  auto make_subconstraint(Vertex_handle va, Vertex_handle vb) {
    return make_sorted_pair(va, vb);
  }

  void add_to_subconstraints_to_conform(Vertex_handle va, Vertex_handle vb,
                                        Constraint_id id) {
    const auto pair = make_subconstraint(va, vb);
#if CGAL_DEBUG_CDT_3 & 32
    std::cerr << "tr.subconstraints_to_conform.push("
              << display_subcstr(pair) << ")\n";
#endif // CGAL_DEBUG_CDT_3
    subconstraints_to_conform.push({pair, id});
  }

  template <typename Visitor>
  Constraint_id insert_constrained_edge_impl(Vertex_handle va, Vertex_handle vb,
                                             Visitor&) {
    if(va != vb) {
      const Constraint_id c_id = constraint_hierarchy.insert_constraint(va, vb);
      // traverse all the vertices along [va, vb] and add pairs of consecutive
      // vertices as sub-constraints.
      std::for_each(tr.segment_traverser_simplices_begin(va, vb), tr.segment_traverser_simplices_end(),
                    [&, prev = Vertex_handle{}](auto simplex) mutable {
                      // std::cerr << "- " << oformat(simplex, With_point_tag{}) << '\n';
                      if(simplex.dimension() == 0) {
                        const auto v = static_cast<Vertex_handle>(simplex);
                        v->c_id = c_id.vl_ptr();
                        ++v->nb_of_incident_constraints;
                        if(prev != Vertex_handle{}) {
                          if(v != vb) {
                            v->set_vertex_type(CDT_3_vertex_type::STEINER_ON_EDGE);
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

  template <typename Visitor>
  Vertex_handle insert_impl(const Point &p, Locate_type lt, Cell_handle c,
                            int li, int lj, Visitor& visitor)
  {
    Vertex_handle v1, v2;
    Vertex_handle new_vertex;
    bool split_constrained_edge = false;
    if(lt == T_3::EDGE) {
      v1 = c->vertex(li);
      v2 = c->vertex(lj);
      if(constraint_hierarchy.is_subconstrained_edge(v1, v2)) {
        split_constrained_edge = true;
      }
    }
    switch (tr.dimension()) {
    case 3: {
      typename T_3::Conflict_tester_3 tester(p, this);
      new_vertex = tr.insert_in_conflict(p, lt, c, li, lj, tester,
                                         visitor);
      break;
    } // dim 3
    case 2: {
      typename T_3::Conflict_tester_2 tester(p, this);
      new_vertex = tr.insert_in_conflict(p, lt, c, li, lj, tester,
                                         visitor);
      break;
    } // dim 2
    default:
      // dimension <= 1
      // Do not use the generic insert.
      new_vertex = tr.insert(p, c);
      break;
    }
    if(split_constrained_edge) {
      constraint_hierarchy.split_constraint(v1, v2, new_vertex);
    }
    return new_vertex;
  }

public:
  Vertex_handle insert(const Point &p, Locate_type lt, Cell_handle c,
                       int li, int lj)
  {
    auto v = insert_impl(p, lt, c, li, lj, insert_in_conflict_visitor);
    restore_Delaunay(insert_in_conflict_visitor);
    return v;
  }

  Vertex_handle insert(const Point &p, Cell_handle start = {}) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr.locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }

  Constraint_id insert_constrained_edge(Vertex_handle va, Vertex_handle vb)
  {
    const auto id = insert_constrained_edge_impl(va, vb, insert_in_conflict_visitor);
    restore_Delaunay(insert_in_conflict_visitor);
    return id;
  }

  bool is_conforming() const {
    return std::all_of(constraint_hierarchy.sc_begin(),
                       constraint_hierarchy.sc_end(),
                       [this](const auto &sc) {
                         const auto va = sc.first.first;
                         const auto vb = sc.first.second;
                         const auto is_edge = this->tr.tds().is_edge(va, vb);
#if CGAL_DEBUG_CDT_3 & 128 && __has_include(<format>)
                         std::cerr << std::format("is_conforming>> Edge is 3D: {}  ({} , {})\n",
                                                  is_edge,
                                                  CGAL::IO::oformat(this->point(va)),
                                                  CGAL::IO::oformat(this->point(vb)));
#endif // CGAL_DEBUG_CDT_3
                         return is_edge;
                       });
  }

  bool write_missing_segments_file(std::ostream &out) {
    bool any_missing_segment = false;
    std::for_each(
        constraint_hierarchy.sc_begin(), constraint_hierarchy.sc_end(),
        [this, &out, &any_missing_segment](const auto &sc) {
          if (!tr.tds().is_edge(sc.first.first, sc.first.second)) {
            const auto v0 = sc.first.first;
            const auto v1 = sc.first.second;
            out << "2 " << this->tr.point(v0) << " " << this->tr.point(v1)
                << '\n';
            any_missing_segment = true;
          }
        });
    return any_missing_segment;
  }

  void write_all_segments_file(std::ostream &out) {
    std::for_each(
        constraint_hierarchy.sc_begin(), constraint_hierarchy.sc_end(),
        [this, &out](const auto &sc) {
          const auto v0 = sc.first.first;
          const auto v1 = sc.first.second;
          out << "2 " << this->tr.point(v0) << " " << this->tr.point(v1) << '\n';
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
  void debug_dump(std::string filename) {
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
    while(!subconstraints_to_conform.empty()) {
      auto pair = subconstraints_to_conform.top();
      subconstraints_to_conform.pop();
      if(!constraint_hierarchy.is_subconstrained_edge(pair.first.first,
                                                      pair.first.second)) {
        continue;
      }
#if CGAL_DEBUG_CDT_3 & 32
      std::cerr << "tr.subconstraints_to_conform.pop()="
                << display_subcstr(pair.first) << "\n";
#endif // CGAL_DEBUG_CDT_3
      conform_subconstraint(pair.first, pair.second, visitor);
    }
  }

  template <typename Visitor>
  Vertex_handle insert_Steiner_point_on_subconstraint(
      Point steiner_pt, Cell_handle hint,
      Subconstraint subconstraint, Constraint_id constraint, Visitor& visitor)
  {
    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    Locate_type lt;
    int li, lj;
    const Cell_handle c = tr.locate(steiner_pt, lt, li, lj, hint);
    const Vertex_handle v = this->insert_impl(steiner_pt, lt, c, li, lj, visitor);
    v->set_vertex_type(CDT_3_vertex_type::STEINER_ON_EDGE);
    if(lt != T_3::VERTEX) {
      v->nb_of_incident_constraints = 1;
      v->c_id = constraint.vl_ptr();
    }
    constraint_hierarchy.add_Steiner(va, vb, v);
    visitor.insert_Steiner_point_on_constraint(constraint, va, vb, v);
    add_to_subconstraints_to_conform(va, v, constraint);
    add_to_subconstraints_to_conform(v, vb, constraint);
    return v;
  }

  /// Return `true` if a Steiner point was inserted
  template <typename Visitor>
  bool conform_subconstraint(Subconstraint subconstraint,
                             Constraint_id constraint,
                             Visitor& visitor)
  {
    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    CGAL_assertion(va != vb);
    if (!tr.tds().is_edge(va, vb)) {
      const auto& [steiner_pt, hint, ref_vertex] = construct_Steiner_point(constraint, subconstraint);
      [[maybe_unused]] const auto v =
          insert_Steiner_point_on_subconstraint(steiner_pt, hint, subconstraint, constraint, visitor);
#ifdef CGAL_DEBUG_CDT_3
      std::cerr << "  new vertex " << display_vert(v) << '\n';
#endif // CGAL_DEBUG_CDT_3
      return true;
    }

    return false;
  }

  Constraint_id constraint_from_extremities(Vertex_handle va, Vertex_handle vb) const {
    if (va->nb_of_incident_constraints == 0 || vb->nb_of_incident_constraints == 0)
    {
      return {};
    }
    Constraint_id c_id = constraint_around(va, vb, false);
    if(c_id != Constraint_id{}) return c_id;
    c_id = constraint_around(vb, va, false);
    if(c_id != Constraint_id{}) return c_id;
    c_id = constraint_around(va, vb, true);
    return c_id;
  }

  auto constraint_extremitites(Constraint_id c_id) const {
      CGAL_assertion(std::find(this->constraint_hierarchy.c_begin(),
                               this->constraint_hierarchy.c_end(), c_id) != this->constraint_hierarchy.c_end());
      CGAL_assertion(this->constraint_hierarchy.vertices_in_constraint_begin(c_id) !=
                     this->constraint_hierarchy.vertices_in_constraint_end(c_id));
#if CGAL_DEBUG_CDT_3 & 8
      std::cerr << "constraint " << (void*) c_id.vl_ptr() << " has "
                << c_id.vl_ptr()->skip_size() << " vertices\n";
#endif // CGAL_DEBUG_CDT_3
      auto it = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
      const auto end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
      const auto c_va = *it;
      while(std::next(it) != end)
        ++it;
      const auto c_vb = *it;
    return std::make_pair(c_va, c_vb);
  }

  Constraint_id constraint_around(Vertex_handle va, Vertex_handle vb, bool expensive = true) const {
    auto constraint_id_goes_to_vb = [this, va, vb](Constraint_id c_id) {
      const auto [c_va, c_vb] = constraint_extremitites(c_id);
      if (va == c_va && vb == c_vb)
        return true;
      if (vb == c_va && va == c_vb)
        return true;
      return false;
    }; // end lambda constraint_id_goes_to_vb
    if (va->nb_of_incident_constraints == 1)
    {
      const Constraint_id c_id = va->constraint_id(*this);
      CGAL_assertion(c_id != Constraint_id{});
      if(constraint_id_goes_to_vb(c_id)) return c_id;
    } else if (expensive == true && va->nb_of_incident_constraints > 1) {
      boost::container::small_vector<Vertex_handle, 64> adj_vertices;
      this->finite_adjacent_vertices(va, std::back_inserter(adj_vertices));
      for(auto other_v: adj_vertices) {
        for(auto context: this->constraint_hierarchy.contexts_range(va, other_v)) {
          const Constraint_id c_id = context.id();
          if(constraint_id_goes_to_vb(c_id)) return c_id;
        }
      }
    }
    return Constraint_id{};
  }

  struct Construct_Steiner_point_return_type {
    typename T_3::Geom_traits::Point_3 point;
    Cell_handle hint;
    Vertex_handle reference_vertex;
  };

  Construct_Steiner_point_return_type
  construct_Steiner_point(Constraint_id constraint_id, Subconstraint subconstraint)
  {
    auto& gt = tr.geom_traits();
    auto angle_functor = gt.angle_3_object();
    auto compare_angle_functor = gt.compare_angle_3_object();
    auto vector_functor = gt.construct_vector_3_object();
    auto midpoint_functor = gt.construct_midpoint_3_object();
    auto scaled_vector_functor = gt.construct_scaled_vector_3_object();
    auto sq_length_functor = gt.compute_squared_length_3_object();
    auto sc_product_functor = gt.compute_scalar_product_3_object();
    auto translate_functor = gt.construct_translated_point_3_object();

    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    const auto& pa = tr.point(va);
    const auto& pb = tr.point(vb);

    if(this->dimension() < 2) {
      std::cerr << "dim < 2: midpoint\n";
      return {midpoint_functor(pa, pb), va->cell(), va};
    }

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "construct_Steiner_point( " << display_vert(va) << " , "
              << display_vert(vb) << " )\n";
#endif // CGAL_DEBUG_CDT_3

#ifdef CGAL_DEBUG_CDT_3
    [[maybe_unused]] auto debug_simplex = [&](auto simplex) {
      std::cerr << " - " << oformat(simplex, With_point_tag{}) << '\n';
    };
#endif // CGAL_DEBUG_CDT_3

    namespace bc = boost::container;
    bc::flat_set<Vertex_handle, std::less<Vertex_handle>,
                 bc::small_vector<Vertex_handle, 256>>
        encroaching_vertices;
    auto register_vertex = [this,&encroaching_vertices](Vertex_handle v) {
      if(tr.is_infinite(v)) return;
      // std::cerr << "register_vertex " << display_vert(v) << '\n';
      encroaching_vertices.insert(v);
    };
    auto fill_encroaching_vertices = [&](const auto simplex) {
#ifdef CGAL_DEBUG_CDT_3
      debug_simplex(simplex);
#endif // CGAL_DEBUG_CDT_3
      auto visit_cell = [&](Cell_handle cell) {
        for(int i = 0, end = this->tr.dimension() + 1; i < end; ++i) {
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
        const auto [other_cell, other_index] = tr.mirror_facet({cell, facet_index});
        register_vertex(other_cell->vertex(other_index));
        break;
      }
      case 1: {
        auto circ = tr.incident_cells(static_cast<Edge>(simplex));
        CGAL_assertion(circ != nullptr);
        const auto end = circ;
        do {
          visit_cell(circ);
        } while(++circ != end);
      } break;
      case 0: {
        const auto v = static_cast<Vertex_handle>(simplex);
        if(v != va && v != vb) {
          std::cerr << "!! The constraint passes through a vertex! ";
          std::cerr << "  -> " << display_vert(v) << '\n';
          debug_dump("bug-through-vertex");
          CGAL_error();
        }
      } break;
      default: CGAL_unreachable();
      } // end switch
    };
    std::for_each(tr.segment_traverser_simplices_begin(va, vb), tr.segment_traverser_simplices_end(),
                  fill_encroaching_vertices);
    auto vector_of_encroaching_vertices = encroaching_vertices.extract_sequence();
#if CGAL_DEBUG_CDT_3 & 16
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
                                                    this->tr.point(v),
                                                    pb) == ACUTE;
                             });
#if CGAL_DEBUG_CDT_3 & 16
    std::cerr << "  -> vector_of_encroaching_vertices (after filter):\n";
    std::for_each(vector_of_encroaching_vertices.begin(), end,
                  [this](Vertex_handle v){
                    std::cerr << "    " << this->display_vert(v) << '\n';
                  });
#endif // CGAL_DEBUG_CDT_3
    CGAL_assertion(vector_of_encroaching_vertices.begin() != end);

    const auto reference_vertex_it = std::max_element(
        vector_of_encroaching_vertices.begin(), end,
        [pa, pb, &compare_angle_functor, this](Vertex_handle v1,
                                               Vertex_handle v2) {
          return compare_angle_functor(pa, this->tr.point(v1), pb, pa,
                                       this->tr.point(v2), pb) == SMALLER;
        });
    CGAL_assertion(reference_vertex_it != end);
#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> reference point: " << display_vert(*reference_vertex_it)
              << '\n';
#endif // CGAL_DEBUG_CDT_3
    const auto reference_vertex = *reference_vertex_it;
    const auto& reference_point = tr.point(reference_vertex);

    if(reference_vertex->is_Steiner_vertex_on_edge()) {
      CGAL_assertion(reference_vertex->nb_of_incident_constraints == 1);
      const auto ref_constraint_id = reference_vertex->constraint_id(*this);
      const auto [ref_va, ref_vb] = constraint_extremitites(ref_constraint_id);
      const auto [orig_va, orig_vb] = constraint_extremitites(constraint_id);
      const auto& orig_pa = tr.point(orig_va);
      const auto& orig_pb = tr.point(orig_vb);
      const auto vector_orig_ab = vector_functor(orig_pa, orig_pb);
      auto return_orig_result_point =
          [&](auto lambda, Point orig_pa, Point orig_pb)
              -> Construct_Steiner_point_return_type
          {
            const auto vector_orig_ab = vector_functor(orig_pa, orig_pb);
            const auto result_point = (lambda < 0.2 || lambda > 0.8)
                                          ? midpoint_functor(pa, pb)
                                          : translate_functor(orig_pa, scaled_vector_functor(vector_orig_ab, lambda));

#ifdef CGAL_DEBUG_CDT_3
            std::cerr << "ref lambda = " << lambda << '\n';
            std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_DEBUG_CDT_3
            return {result_point, reference_vertex->cell(), reference_vertex};
          };

      const auto length_orig_ab = CGAL::sqrt(sq_length_functor(vector_orig_ab));
      if(ref_va == orig_va || ref_vb == orig_va) {
        const auto vector_orig_a_ref = vector_functor(orig_pa, reference_point);
        const auto length_orig_a_ref = CGAL::sqrt(sq_length_functor(vector_orig_a_ref));
        const auto lambda = length_orig_a_ref / length_orig_ab;
        return return_orig_result_point(lambda, orig_pa, orig_pb);
      } else if(ref_va == orig_vb || ref_vb == orig_vb) {
        const auto vector_orig_b_ref = vector_functor(orig_pb, reference_point);
        const auto length_orig_b_ref = CGAL::sqrt(sq_length_functor(vector_orig_b_ref));
        const auto lambda = length_orig_b_ref / length_orig_ab;
        return return_orig_result_point(lambda, orig_pb, orig_pa);
      }
    }
    // compute the projection of the reference point
    const auto vector_ab = vector_functor(pa, pb);
    const auto vector_a_ref = vector_functor(pa, reference_point);
    const auto lambda = sc_product_functor(vector_a_ref, vector_ab) / sq_length_functor(vector_ab);
    const auto result_point = (lambda < 0.2 || lambda > 0.8)
                                  ? midpoint_functor(pa, pb)
                                  : translate_functor(pa, scaled_vector_functor(vector_ab, lambda));

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "lambda = " << lambda << '\n';
    std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_DEBUG_CDT_3
    return {result_point, reference_vertex->cell(), reference_vertex};
  }

protected:
  T_3& tr = *this;
  Compare_vertex_handle comp = {this};
  Constraint_hierarchy constraint_hierarchy = {comp};
  Insert_in_conflict_visitor insert_in_conflict_visitor = {*this};

  std::stack<std::pair<Subconstraint, Constraint_id> >
    subconstraints_to_conform;
};

} // end CGAL

#endif // CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
