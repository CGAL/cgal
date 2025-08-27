#ifndef HEXMESHING_LAPLACIAN_SMOOTHING_H
#define HEXMESHING_LAPLACIAN_SMOOTHING_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_mark_utils.h>
#include <CGAL/Eigen_matrix.h>
#include <vector>

namespace CGAL::internal::Hexmeshing {
  // laplacian smoothing on points which isn't on the surface of the grid and don't have mark
  void laplacian_smoothing_for_unmarked_cells(LCC& lcc, size_type surface_mark) {
    auto vertices = lcc.one_dart_per_cell<0>();
    auto edges = lcc.one_dart_per_cell<1>();
    auto faces = lcc.one_dart_per_cell<2>();
    size_type side_mark = lcc.get_new_mark();

    const int count_vertices = set_vertex_ids(lcc);
    std::vector<Vector> P_new(count_vertices, CGAL::NULL_VECTOR);
    std::vector<int> count(count_vertices);
    for(auto edge = edges.begin(); edge != edges.end(); edge++) {
      auto redge = lcc.beta<1>(edge);
      int id1 = lcc.attribute<0>(edge)->id, id2 = lcc.attribute<0>(redge)->id;
      P_new[id1] += (lcc.point(redge) - CGAL::ORIGIN);
      P_new[id2] += (lcc.point(edge) - CGAL::ORIGIN);
      count[id1]++;
      count[id2]++;
    }
    for(auto face = faces.begin(); face != faces.end(); face++) {
      if(lcc.is_free<3>(face))
        mark_k_cells_of_i_cell<2, 0>(lcc, face, side_mark);
    }

    for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
      if(lcc.is_marked(vertex, surface_mark) or lcc.is_marked(vertex, side_mark)) continue;

      int id = lcc.attribute<0>(vertex)->id;

      assert(count[id]);

      lcc.point(vertex) = CGAL::ORIGIN + P_new[id]/count[id];
    }

    lcc.free_mark(side_mark);
  }

  void volume_smoothing(LCC& lcc, size_type surface_mark) {
    for(int _ = 2; _--;)
      laplacian_smoothing_for_unmarked_cells(lcc, surface_mark);
    // for(int _ = 2; _--;)
    //   shape_improvement_for_unmarked_cells(lcc, surface_mark);
  }

  std::vector<std::vector<Dart_handle>> get_neighbors_list_for_smoothing(LCC& lcc, size_type surface_mark, size_type inner_mark) {
    int count_vertices = set_vertex_ids(lcc);
    std::vector<std::vector<Dart_handle>> neighbors_list(count_vertices);

    auto edges = lcc.one_dart_per_cell<1>();
    for(auto edge = edges.begin(); edge != edges.end(); edge++) {
      bool inside = lcc.is_marked(edge, inner_mark);
      bool outside = !inside;
      Dart_handle e = lcc.beta(edge, 3, 2);

      while(e != nullptr and e != edge and lcc.attribute<3>(e) != nullptr) {
        if(lcc.is_marked(e, inner_mark)) {
          inside = true;
        }
        else {
          outside = true;
        }
        e = lcc.beta(e, 3, 2);
      }
      if(e == nullptr or lcc.attribute<3>(e) == nullptr) {
        outside = true;
        e = lcc.beta(edge, 2);
        while(e != nullptr and lcc.attribute<3>(e) != nullptr) {
          if(lcc.is_marked(e, inner_mark)) {
            inside = true;
          }
          e = lcc.beta(e, 3, 2);
        }
      }
      if(inside and outside and lcc.is_whole_cell_marked<1>(edge, surface_mark)) {
        int id1 = lcc.attribute<0>(edge)->id;
        int id2 = lcc.attribute<0>(lcc.beta<1>(edge))->id;
        neighbors_list[id1].emplace_back(lcc.beta<1>(edge));
        neighbors_list[id2].emplace_back(edge);
      }
    }

    return neighbors_list;
  }

  // normals for vertices need to be set
  void surface_smoothing(LCC& lcc, size_type surface_mark, size_type inner_mark, const double ridge_ratio=0.005) {
    auto vertices = lcc.one_dart_per_cell<0>();

    std::vector<std::vector<Dart_handle>> neighbors_list = get_neighbors_list_for_smoothing(lcc, surface_mark, inner_mark);

    std::vector<std::pair<Dart_handle, Point>> new_points;

    for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
      auto attr = lcc.attribute<0>(vertex);
      int id = attr->id;
      if(neighbors_list[id].empty()) continue;

      Vector N_k = attr->normal;
      auto &&[T_1, T_2] = get_orthogonal_vectors(N_k);

      CGAL::Eigen_matrix<double, 5, 5> mat;
      CGAL::Eigen_vector<double, 5> vec;
      for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
          mat.matrix()(i, j) = 0.0;
        }
        vec(i) = 0.0;
      }
      Vector barycenter = CGAL::NULL_VECTOR;
      
      for(auto neighbor: neighbors_list[id]) {
        Vector diff = lcc.point(neighbor) - lcc.point(vertex);
        double x = diff * T_1, y = diff * T_2, z = diff * N_k;
        std::array<double, 5> qd = {x, y, x*x, x*y, y*y};
        double w = 1. / std::sqrt(diff * diff);
        for(int i = 0; i < 5; i++) {
          for(int j = 0; j < 5; j++) {
            mat.matrix()(i, j) += w*qd[i]*qd[j];
          }
          vec(i) = vec(i) + w*z*qd[i];
        }
        barycenter += Vector(x, y, z);
      }
      barycenter /= neighbors_list[id].size();
      
      double lam = 0.;
      for(int i = 0; i < 5; i++) lam += mat(i, i);
      // the main diagonal values supposed here are mat(2, 2) and mat(4, 4)
      // I recommend ridge_ratio to be 1/200 so that lambda becomes 1/100 of the main values
      lam *= ridge_ratio;
      for(int i = 0; i < 5; i++) mat.matrix()(i, i) += lam;
      auto a_k = mat.ldlt().solve(vec);
      auto Q_k = [&](double x, double y) -> double {
        return a_k[0]*x + a_k[1]*y + a_k[2]*x*x + a_k[3]*x*y + a_k[4]*y*y;
      };
      auto new_point = lcc.point(vertex) + barycenter.x()*T_1 + barycenter.y()*T_2 + Q_k(barycenter.x(), barycenter.y())*N_k;
      new_points.emplace_back(vertex, new_point);
    }

    for(auto &[dart, new_point]: new_points) {
      lcc.point(dart) = new_point;
    }
  }
}


#endif