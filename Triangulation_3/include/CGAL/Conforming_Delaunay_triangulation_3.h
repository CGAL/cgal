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

  static std::string io_signature() {
    return Get_io_signature<Vb>()();
  }
};


template <typename T_3>
class Conforming_Delaunay_triangulation_3 : public T_3 {
public:
  using Vertex_handle = typename T_3::Vertex_handle;
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
    os << oformat(v) << "=(" << this->tr.point(v) << ")";
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
  };

  auto make_subconstraint(Vertex_handle va, Vertex_handle vb) {
    return make_sorted_pair(va, vb);
  }

  void add_to_subconstraints_to_conform(Vertex_handle va, Vertex_handle vb,
                                        Constraint_id id) {
    const auto pair = make_subconstraint(va, vb);
#if CGAL_DEBUG_CDT_3
    std::cerr << "tr.subconstraints_to_conform.push("
              << display_subcstr(pair) << ")\n";
#endif // CGAL_DEBUG_CDT_3
    subconstraints_to_conform.push({pair, id});
  }

  template <typename Visitor>
  Constraint_id insert_constrained_edge_impl(Vertex_handle va, Vertex_handle vb,
                                             Visitor& visitor) {
    if(va != vb) {
      const Constraint_id c_id = constraint_hierarchy.insert_constraint(va, vb);
      va->c_id = c_id.vl_ptr();
      vb->c_id = c_id.vl_ptr();
      ++va->nb_of_incident_constraints;
      ++vb->nb_of_incident_constraints;
      add_to_subconstraints_to_conform(va, vb, c_id);
      restore_Delaunay(visitor);
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
      debug_dump("dump-bug");
      CGAL_error_msg("case not yet handled");
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
    return insert_constrained_edge_impl(va, vb, insert_in_conflict_visitor);
  }

  bool is_conforming() const {
    return std::all_of(constraint_hierarchy.sc_begin(),
                       constraint_hierarchy.sc_end(),
                       [this](const auto &sc) {
                         return this->tr.tds().is_edge(sc.first.first,
                                                       sc.first.second);
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
      CGAL::Mesh_3::save_binary_file(dump, *this);
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
#if CGAL_DEBUG_CDT_3
      std::cerr << "tr.subconstraints_to_conform.pop()="
                << display_subcstr(pair.first) << "\n";
#endif // CGAL_DEBUG_CDT_3
      conform_subconstraint(pair.first, pair.second, visitor);
    }
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
      const auto& [steiner_pt, hint] = construct_Steiner_point(subconstraint);
      Locate_type lt;
      int li, lj;
      const Cell_handle c = tr.locate(steiner_pt, lt, li, lj, hint);
      const Vertex_handle v = this->insert_impl(steiner_pt, lt, c, li, lj,
                                                visitor);
      if(lt != T_3::VERTEX) {
        v->nb_of_incident_constraints = 1;
        v->c_id = constraint.vl_ptr();
#if CGAL_DEBUG_CDT_3
        std::cerr << "New Steiner vertex (#" << tr.number_of_vertices() << "): "
                  << display_vert(v) << '\n';
#endif // CGAL_DEBUG_CDT_3
      }
      CGAL_assertion(v != va);
      CGAL_assertion(v != vb);
      constraint_hierarchy.add_Steiner(va, vb, v);
      add_to_subconstraints_to_conform(va, v, constraint);
      add_to_subconstraints_to_conform(v, vb, constraint);
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

  Constraint_id constraint_around(Vertex_handle va, Vertex_handle vb, bool expensive = true) const {
    auto constraint_id_goes_to_vb = [this, va, vb](Constraint_id c_id) {
      CGAL_assertion(std::find(this->constraint_hierarchy.c_begin(),
                               this->constraint_hierarchy.c_end(), c_id) != this->constraint_hierarchy.c_end());
      CGAL_assertion(this->constraint_hierarchy.vertices_in_constraint_begin(c_id) !=
                     this->constraint_hierarchy.vertices_in_constraint_end(c_id));
#ifdef CGAL_DEBUG_CDT_3
      std::cerr << "constraint " << (void*) c_id.vl_ptr() << " has "
                << c_id.vl_ptr()->skip_size() << " vertices\n";
#endif
      auto it = this->constraint_hierarchy.vertices_in_constraint_begin(c_id);
      const auto end = this->constraint_hierarchy.vertices_in_constraint_end(c_id);
      const auto c_va = *it;
      while(std::next(it) != end)
        ++it;
      const auto c_vb = *it;
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
  };

  auto construct_Steiner_point(Subconstraint subconstraint) -> Construct_Steiner_point_return_type
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
      return {midpoint_functor(pa, pb), va->cell()};
    }

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "construct_Steiner_point( " << display_vert(va) << " , "
              << display_vert(vb) << " )\n";
#endif // CGAL_DEBUG_CDT_3
    const CGAL::Triangulation_segment_cell_iterator_3<T_3> cell_traverser_begin{&tr, va, vb};
    const auto cell_traverser_end = cell_traverser_begin.end();

    namespace bc = boost::container;
    bc::flat_set<Vertex_handle, std::less<Vertex_handle>,
                 bc::small_vector<Vertex_handle, 256>>
        encroaching_vertices;
    auto fill_encroaching_vertices = [this,
                                      &encroaching_vertices](const auto &cell) {
      for (int i = 0, end = this->tr.dimension(); i < end; ++i) {
        const auto v = cell.vertex(i);
        encroaching_vertices.insert(v);
      }
    };
    std::for_each(cell_traverser_begin, cell_traverser_end,
                  fill_encroaching_vertices);
    auto vector_of_encroaching_vertices = encroaching_vertices.extract_sequence();
#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> vector_of_encroaching_vertices (before filter):\n";
    std::for_each(vector_of_encroaching_vertices.begin(),
                  vector_of_encroaching_vertices.end(),
                  [this](Vertex_handle v){
                    std::cerr << "    " << this->display_vert(v) << '\n';
                  });
#endif
    auto end = std::remove_if(vector_of_encroaching_vertices.begin(),
                              vector_of_encroaching_vertices.end(),
                              [pa, pb, &angle_functor, this](Vertex_handle v) {
                               return angle_functor(pa,
                                                    this->tr.point(v),
                                                    pb) == ACUTE;
                             });
#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> vector_of_encroaching_vertices (after filter):\n";
    std::for_each(vector_of_encroaching_vertices.begin(), end,
                  [this](Vertex_handle v){
                    std::cerr << "    " << this->display_vert(v) << '\n';
                  });
#endif
    CGAL_assertion(vector_of_encroaching_vertices.begin() != end);

    auto reference_point_it = std::max_element(
        vector_of_encroaching_vertices.begin(), end,
        [pa, pb, &compare_angle_functor, this](Vertex_handle v1,
                                               Vertex_handle v2) {
          return compare_angle_functor(pa, this->tr.point(v1), pb, pa,
                                       this->tr.point(v2), pb) == SMALLER;
        });
    CGAL_assertion(reference_point_it != end);
#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> reference point: " << display_vert(*reference_point_it)
              << '\n';
#endif // CGAL_DEBUG_CDT_3
    const auto &reference_point = tr.point(*reference_point_it);
    const auto cell_incident_to_reference_point = (*reference_point_it)->cell();

    // compute the projection of the reference point
    const auto vector_ab = vector_functor(pa, pb);
    const auto vector_a_ref = vector_functor(pa, reference_point);
    const auto lambda = sc_product_functor(vector_a_ref, vector_ab) /
                        CGAL::sqrt(sq_length_functor(vector_ab));

    const auto result_point =
     (lambda < 0.2 || lambda > 0.8) ?
      midpoint_functor(pa, pb) :
      translate_functor(pa, scaled_vector_functor(vector_ab, lambda));

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_DEBUG_CDT_3
    return { result_point, cell_incident_to_reference_point };
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
