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

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp> /// @TODO Requires Boost 1.66

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

};


template <typename T_3>
class Triangulation_conformer_3 : public T_3 {
  using Vertex_handle = typename T_3::Vertex_handle;
  using Cell_handle = typename T_3::Cell_handle;
  using Point = typename T_3::Point;
  using Locate_type = typename T_3::Locate_type;

  struct Compare_vertex_handle {
    const T_3* tr;
    Compare_vertex_handle(const T_3* tr) : tr(tr) {}
    bool operator()(const Vertex_handle va, const Vertex_handle vb) const {
      return tr->compare_xyz(tr->point(va), tr->point(vb)) == SMALLER;
    }
  };

  using Constraint_hierarchy =
      Polyline_constraint_hierarchy_2<Vertex_handle, Compare_vertex_handle, const Point&>;
  using Constraint_id = typename Constraint_hierarchy::Constraint_id;
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
    Triangulation_conformer_3<T_3>& self;
  public:
    Insert_in_conflict_visitor(Triangulation_conformer_3& self) : self(self) {}

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

  Vertex_handle private_insert(const Point &p, Locate_type lt, Cell_handle c,
                               int li, int lj)
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
                                         insert_in_conflict_visitor);
      break;
    } // dim 3
    case 2: {
      typename T_3::Conflict_tester_2 tester(p, this);
      new_vertex = tr.insert_in_conflict(p, lt, c, li, lj, tester,
                                         insert_in_conflict_visitor);
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
      CGAL_error_msg("case not yet handled");
    }
    return new_vertex;
  }

public:
  Triangulation_conformer_3()
      : tr(*this), comp(&tr), constraint_hierarchy(comp) {
  }

  Vertex_handle insert(const Point &p, Locate_type lt, Cell_handle c,
                               int li, int lj)
  {
    auto v = private_insert(p, lt, c, li, lj);
    restore_Delaunay();
    return v;
  }

  Vertex_handle insert(const Point &p, Cell_handle start = {}) {
    Locate_type lt;
    int li, lj;

    Cell_handle c = tr.locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }

  void insert_constrained_edge(Vertex_handle va, Vertex_handle vb) {
    if(va != vb) {
      const Constraint_id c_id = constraint_hierarchy.insert_constraint(va, vb);
      va->c_id = c_id.vl_ptr();
      vb->c_id = c_id.vl_ptr();
      ++va->nb_of_incident_constraints;
      ++vb->nb_of_incident_constraints;
      add_to_subconstraints_to_conform(va, vb, c_id);
      restore_Delaunay();
    }
  }

  bool is_conforming() const {
    return std::all_of(constraint_hierarchy.sc_begin(),
                       constraint_hierarchy.sc_end(),
                       [this](const auto &sc) {
                         return this->tr.tds().is_edge(sc.first.first,
                                                       sc.first.second);
                       });
  }

protected:
  void restore_Delaunay() {
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
      conform_subconstraint(pair.first, pair.second);
    }
  }

  /// Return `true` if a Steiner point was inserted
  bool conform_subconstraint(Subconstraint subconstraint,
                             Constraint_id constraint)
  {
    const Vertex_handle va = subconstraint.first;
    const Vertex_handle vb = subconstraint.second;
    CGAL_triangulation_assertion(va != vb);
    if (!tr.tds().is_edge(va, vb)) {
      const auto& [steiner_pt, hint] = construct_Steiner_point(subconstraint);
      Locate_type lt;
      int li, lj;
      const Cell_handle c = tr.locate(steiner_pt, lt, li, lj, hint);
      const Vertex_handle v = this->private_insert(steiner_pt, lt, c, li, lj);
      if(lt != T_3::VERTEX) {
        v->nb_of_incident_constraints = 1;
        v->c_id = constraint.vl_ptr();
#if CGAL_DEBUG_CDT_3
        std::cerr << "New Steiner vertex: " << display_vert(v) << '\n';
#endif // CGAL_DEBUG_CDT_3
      }
      CGAL_triangulation_assertion(v != va);
      CGAL_triangulation_assertion(v != vb);
      constraint_hierarchy.add_Steiner(va, vb, v);
      add_to_subconstraints_to_conform(va, v, constraint);
      add_to_subconstraints_to_conform(v, vb, constraint);
      return true;
    }

    return false;
  }

  auto construct_Steiner_point(Subconstraint subconstraint)
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

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "construct_Steiner_point( " << display_vert(va) << " , "
              << display_vert(vb) << " )\n";
#endif // CGAL_DEBUG_CDT_3
    const CGAL::Triangulation_segment_cell_iterator_3<T_3> cell_traverser_begin{tr, va, vb};
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
                  [=](Vertex_handle v){
                    std::cerr << "    " << display_vert(v) << '\n';
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
                  [=](Vertex_handle v){
                    std::cerr << "    " << display_vert(v) << '\n';
                  });
#endif
    CGAL_triangulation_assertion(vector_of_encroaching_vertices.begin() != end);

    auto reference_point_it = std::max_element(
        vector_of_encroaching_vertices.begin(), end,
        [pa, pb, &compare_angle_functor, this](Vertex_handle v1,
                                               Vertex_handle v2) {
          return compare_angle_functor(pa, this->tr.point(v1), pb, pa,
                                       this->tr.point(v2), pb) == SMALLER;
        });
    CGAL_triangulation_assertion(reference_point_it != end);
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
                        sq_length_functor(vector_ab);

    const auto result_point =
     (lambda < 0.2 || lambda > 0.8) ?
      midpoint_functor(pa, pb) :
      translate_functor(pa, scaled_vector_functor(vector_ab, lambda));

#ifdef CGAL_DEBUG_CDT_3
    std::cerr << "  -> Steiner point: " << result_point << '\n';
#endif // CGAL_DEBUG_CDT_3
    return std::make_pair(result_point, cell_incident_to_reference_point);
  }

private:
  T_3& tr;
  Compare_vertex_handle comp;
  Constraint_hierarchy constraint_hierarchy;
  Insert_in_conflict_visitor insert_in_conflict_visitor = *this;

  std::stack<std::pair<Subconstraint, Constraint_id> >
    subconstraints_to_conform;
};
}

#endif // CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
