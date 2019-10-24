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

#include <CGAL/Triangulation_2/internal/Constraint_hierarchy_2.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp> /// @TODO Requires Boost 1.66

    template <typename T> struct Display_type;

namespace CGAL {

template <typename T_3>
class Triangulation_conformer_3 {
  using Vertex_handle = typename T_3::Vertex_handle;
  using Cell_handle = typename T_3::Cell_handle;
  using Point = typename T_3::Point;

  struct Compare_vertex_handle {
    const T_3* tr;
    Compare_vertex_handle(const T_3* tr) : tr(tr) {}
    bool operator()(const Vertex_handle va, const Vertex_handle vb) const {
      return tr->compare_xyz(tr->point(va), tr->point(vb)) == SMALLER;
    }
  };

  using Constraints_hierarchy =
      Constraint_hierarchy_2<Vertex_handle, Compare_vertex_handle, bool>;

public:
  Triangulation_conformer_3(T_3 &tr)
      : tr(tr), comp(&tr), constraints_hierarchy(comp) {
  }

  void insert_constrained_edge(Vertex_handle va, Vertex_handle vb) {
    if(va != vb)
      constraints_hierarchy.insert_constraint(va, vb);
  }

  void restore_Delaunay() {
    bool not_yet_restored = false;
    do {
      not_yet_restored = false;
      for (auto h_sc : make_range(constraints_hierarchy.sc_begin(),
                                  constraints_hierarchy.sc_end())) {
        const Vertex_handle va = h_sc.first.first;
        const Vertex_handle vb = h_sc.first.second;
        CGAL_triangulation_assertion(va != vb);
        Cell_handle c;
        int i;
        int j;
        if (!tr.is_edge(va, vb, c, i, j)) {
          const auto& [steiner_pt, hint] = construct_Steiner_point(h_sc);
          const Vertex_handle v = tr.insert(steiner_pt, hint);
          CGAL_triangulation_assertion(v != va);
          CGAL_triangulation_assertion(v != vb);
          constraints_hierarchy.add_Steiner(va, vb, v);
          not_yet_restored = true;
          break;
        }
      }
    }
    while(not_yet_restored);
  }
protected:
  auto construct_Steiner_point(typename Constraints_hierarchy::H_sc_to_c_map::value_type h_sc)
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

    const Vertex_handle va = h_sc.first.first;
    const Vertex_handle vb = h_sc.first.second;
    const auto& pa = tr.point(va);
    const auto& pb = tr.point(vb);

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
    auto end = std::remove_if(vector_of_encroaching_vertices.begin(),
                             vector_of_encroaching_vertices.end(),
                              [pa, pb, &angle_functor, this](Vertex_handle v) {
                               return angle_functor(pa,
                                                    this->tr.point(v),
                                                    pb) == ACUTE;
                             });
    CGAL_triangulation_assertion(vector_of_encroaching_vertices.begin() != end);

    auto reference_point_it = std::max_element(
        vector_of_encroaching_vertices.begin(), end,
        [pa, pb, &compare_angle_functor, this](Vertex_handle v1,
                                               Vertex_handle v2) {
          return compare_angle_functor(pa, this->tr.point(v1), pb, pa,
                                       this->tr.point(v2), pb) == SMALLER;
        });
    CGAL_triangulation_assertion(reference_point_it != end);
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

    return std::make_pair(result_point, cell_incident_to_reference_point);
  }

private:
  T_3& tr;
  Compare_vertex_handle comp;
  Constraints_hierarchy constraints_hierarchy;
};

}

#endif // CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
