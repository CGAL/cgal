// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>
//
#ifndef HEXMESHING_MOVE_POINTS_ONTO_MESH_H
#define HEXMESHING_MOVE_POINTS_ONTO_MESH_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_mark_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_resolve_non_manifold_case.h>
#include <CGAL/hexmeshing/Hexmeshing_function_generator.h>
#include <vector>

namespace CGAL::internal::Hexmeshing {
  // func sets intersection and normal for dual edge
  // normal should be normalized
  void move_points_onto_mesh(LCC& lcc, size_type move_mark, DetectingFunction func) {
    set_dual_edges(lcc);

    int count_vertices = set_vertex_ids(lcc);

    auto faces = lcc.one_dart_per_cell<2>();
    std::vector<Vector> P_news(count_vertices, CGAL::NULL_VECTOR), N_news(count_vertices, CGAL::NULL_VECTOR);
    std::vector<int> count_intersect(count_vertices, 0);
    for(auto it = faces.begin(); it != faces.end(); it++) {
      if(func(lcc, it)) {
        auto &face_attr = lcc.attribute<2>(it)->info();
        for(auto vertex = lcc.one_dart_per_incident_cell<0, 2>(it).begin(),
                  end = lcc.one_dart_per_incident_cell<0, 2>(it).end();
                vertex != end; vertex++) {
          int id = lcc.attribute<0>(vertex)->id;
          Vector P_new_add = ((lcc.point(vertex) - CGAL::ORIGIN) - (face_attr.normal * (lcc.point(vertex) - face_attr.intersection))*face_attr.normal);
          P_news[id] += P_new_add;
          
          if(count_intersect[id]) {
            N_news[id] += face_attr.normal;
          }
          else {
            N_news[id] = face_attr.normal;
          }

          count_intersect[id]++;
        }
        mark_k_cells_of_i_cell<2, 0>(lcc, it, move_mark);
      }
    }

    auto vertices = lcc.one_dart_per_cell<0>();
    for(auto it = vertices.begin(); it != vertices.end(); it++) {
      if(lcc.is_marked(it, move_mark)) {
        int id = lcc.attribute<0>(it)->id;
        P_news[id] /= count_intersect[id];
        lcc.point(it) = CGAL::ORIGIN + P_news[id];

        N_news[id] /= count_intersect[id];
        lcc.attribute<0>(it)->normal = N_news[id];
      }
    }
  }

  // void move_points_onto_mesh_with_volume_fraction(LCC& lcc, size_type move_mark, size_type inner_mark, double length_of_4_template, MarkingFunction cellIdentifier, DecideInsideFunction decideFunc) {
  //   set_fraction(lcc, length_of_4_template, cellIdentifier, decideFunc);
  void move_points_onto_mesh_with_volume_fraction(LCC& lcc, size_type move_mark, size_type inner_mark) {
    const double s = 0.5;

    auto volumes = lcc.one_dart_per_cell<3>();
    int cnt = 0;
    for(auto it = volumes.begin(); it != volumes.end(); it++) {
      if(lcc.attribute<3>(it)->info().fraction > s) lcc.mark_cell<3>(it, inner_mark), cnt++;
    }
    // std::cout << "Number of inner volumes is: " << cnt << std::endl;
    resolve_non_manifold_case(lcc, s, inner_mark);

    // normal を全ての dual node にセットして、それを使って p.7 式 (2, 3) を detectIntersection で実装していく

    size_type set_gradient_mark = lcc.get_new_mark();

    auto detectIntersection = detect_intersection_with_volume_fraction(s, inner_mark, set_gradient_mark);
    move_points_onto_mesh(lcc, move_mark, detectIntersection);

    lcc.free_mark(set_gradient_mark);
  }
}



#endif