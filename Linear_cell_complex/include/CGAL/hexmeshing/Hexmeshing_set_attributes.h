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
#ifndef HEXMESHING_SET_ATTRIBUTES_H
#define HEXMESHING_SET_ATTRIBUTES_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>

namespace CGAL::internal::Hexmeshing {
  void __set_centroid(LCC& lcc, Dart_handle dart) {
    auto attr = get_or_create_attr<3>(lcc, dart);
    auto &vol_attr = attr->info();

    Vector centroid = CGAL::NULL_VECTOR;
    int vertex_count = 0;
    
    for(auto it = lcc.one_dart_per_incident_cell<0, 3>(dart).begin(),
              end = lcc.one_dart_per_incident_cell<0, 3>(dart).end(); 
            it != end; it++) {
      centroid = centroid + (lcc.point(it) - CGAL::ORIGIN);
      vertex_count++;
    }

    assert(vertex_count == 8);
    
    vol_attr.centroid = CGAL::ORIGIN + centroid / vertex_count;
  }

  void set_centroids(LCC& lcc) {
    auto volumes = lcc.one_dart_per_cell<3>();
    for(auto it = volumes.begin(); it != volumes.end(); it++) {
      __set_centroid(lcc, it);
    }
  }

  void __set_dual_edge(LCC& lcc, Dart_handle dart) {
    auto attr = get_or_create_attr<2>(lcc, dart);
    auto &face_attr = attr->info();

    Point p1 = lcc.attribute<3>(dart)->info().centroid, p2(0, 0, 0);
    if(lcc.is_free<3>(dart)) {
      int vertex_count = 0;
      for(auto it = lcc.one_dart_per_incident_cell<0, 2>(dart).begin(),
                end = lcc.one_dart_per_incident_cell<0, 2>(dart).end();
              it != end; it++) {
        p2 = p2 + (lcc.point(it) - CGAL::ORIGIN);
        vertex_count++;
      }

      assert(vertex_count == 4);

      p2 = CGAL::ORIGIN + (p2 - CGAL::ORIGIN) / vertex_count;
    }
    else {
      p2 = lcc.attribute<3>(lcc.beta<3>(dart))->info().centroid;
    }

    face_attr.dual_edge = {p1, p2};
  }

  void set_dual_edges(LCC& lcc) {
    set_centroids(lcc);
    auto faces = lcc.one_dart_per_cell<2>();
    for(auto it = faces.begin(); it != faces.end(); it++) {
      __set_dual_edge(lcc, it);
    }
  }

  int set_vertex_ids(LCC& lcc) {
    auto vertices = lcc.one_dart_per_cell<0>();
    int count_vertices = 0;
    for(auto it = vertices.begin(); it != vertices.end(); it++) {
      get_or_create_attr<0>(lcc, it)->id = count_vertices++;
    }
    return count_vertices;
  }

  void __set_fraction(LCC& lcc, Dart_handle dart, int number_of_random_points, RandomPointGenerator& gen, MarkingFunction cellIdentifier, DecideInsideFunction decideFunc) {
    auto attr = get_or_create_attr<3>(lcc, dart);
    auto &vol_attr = attr->info();
    __set_centroid(lcc, dart);

    if(!cellIdentifier(lcc, dart)) {
      vol_attr.fraction = decideFunc(vol_attr.centroid) ? 1.0 : 0.0;
    }
    else {
      std::vector<Point> random_points;
      copy_n(gen, number_of_random_points, std::back_inserter(random_points));
      int count_inner_points = 0;
      for(Point& p: random_points) {
        if(decideFunc(vol_attr.centroid + (p - CGAL::ORIGIN))) count_inner_points++;
      }
      vol_attr.fraction = static_cast<double>(count_inner_points) / number_of_random_points;
    }
  }

  template<int numberOfRandomPoints=101>
  void set_fraction(LCC& lcc, double length_of_4_template, MarkingFunction cellIdentifier, DecideInsideFunction decideFunc) {
    auto volumes = lcc.one_dart_per_cell<3>();
    RandomPointGenerator gen(length_of_4_template * 0.5);
    for(auto it = volumes.begin(); it != volumes.end(); it++) {
      __set_fraction(lcc, it, numberOfRandomPoints, gen, cellIdentifier, decideFunc);
    }
  }

  // 2.2 Estimate gradients at cell centers 参照で作成
  // centroid is assumed to be set
  void __set_gradient_at_dual_node(LCC& lcc, Dart_handle dart) {
    auto neighbors = cells_26_connectivity(lcc, dart);
    CGAL::Eigen_matrix<double, 3, 3> mat;
    CGAL::Eigen_vector<double, 3> vec;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        mat.matrix()(i, j) = 0.0;
      }
      vec(i) = 0.0;
    }

    for(auto neighbor: neighbors) {
      Vector diff = lcc.attribute<3>(neighbor)->info().centroid - lcc.attribute<3>(dart)->info().centroid;
      std::array<double, 3> qd = {diff.x(), diff.y(), diff.z()};
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          mat.matrix()(i, j) += qd[i]*qd[j];
        }
        vec(i) += (lcc.attribute<3>(neighbor)->info().fraction - lcc.attribute<3>(dart)->info().fraction) * qd[i];
      }
    }

    auto grad_frac = mat.ldlt().solve(vec);
    lcc.attribute<3>(dart)->info().gradient = Vector(grad_frac[0], grad_frac[1], grad_frac[2]);
  }
}






#endif