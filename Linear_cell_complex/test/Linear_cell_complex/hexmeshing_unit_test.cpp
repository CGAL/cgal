#include <CGAL/hexmeshing/Hexmeshing_set_attributes.h>
#include <CGAL/hexmeshing/Hexmeshing_laplacian_smoothing.h>
#include <CGAL/hexmeshing/Hexmeshing_move_points_onto_mesh.h>
#include <CGAL/hexmeshing/Hexmeshing_resolve_non_manifold_case.h>

#include "lcc_jacobian.h"
#include "unit_test_framework.h"

#include <vector>
#include <array>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <fstream>



template<int SKBT>
bool test_set_dual_edges();

template<>
bool test_set_dual_edges<1>() {
  LCC lcc = create_grid_mesh(1, 1, 2);
  set_dual_edges(lcc);
  auto faces = lcc.one_dart_per_cell<2>();
  std::vector<Segment> result;
  for(auto face = faces.begin(); face != faces.end(); face++) {
    result.emplace_back(lcc.attribute<2>(face)->info().dual_edge);
  }

  const std::array<Segment, 10> answer_surface = {{Segment(Point(0.5, 0.5, 0.5), Point(0, 0.5, 0.5)), Segment(Point(0.5, 0.5, 0.5), Point(0.5, 0, 0.5)), Segment(Point(0.5, 0.5, 0.5), Point(0.5, 0.5, 0)), Segment(Point(0.5, 0.5, 0.5), Point(1, 0.5, 0.5)), Segment(Point(0.5, 0.5, 0.5), Point(0.5, 1, 0.5)),
                                        Segment(Point(0.5, 0.5, 1.5), Point(0, 0.5, 1.5)), Segment(Point(0.5, 0.5, 1.5), Point(0.5, 0, 1.5)), Segment(Point(0.5, 0.5, 1.5), Point(1, 0.5, 1.5)), Segment(Point(0.5, 0.5, 1.5), Point(0.5, 1, 1.5)), Segment(Point(0.5, 0.5, 1.5), Point(0.5, 0.5, 2))}};
  const std::array<Segment, 1> answer_inner = {{Segment(Point(0.5, 0.5, 0.5), Point(0.5, 0.5, 1.5))}};

  if(result.size() != answer_surface.size() + answer_inner.size())
    return false;

  for(auto answer: answer_surface) {
    Point s = answer.source();
    Point t = answer.target();
    bool exist = false;
    for(auto seg: result) {
      Point s_ = seg.source();
      Point t_ = seg.target();
      if(is_approximately_equal_Point(s, s_) and is_approximately_equal_Point(t, t_)) {
        exist = true;
        break;
      }
    }
    if(!exist)
      return false;
  }
  
  for(auto answer: answer_inner) {
    bool exist = false;
    for(auto seg: result) {
      if(is_approximately_equal_Segment(answer, seg)) {
        exist = true;
        break;
      }
    }
    if(!exist)
      return false;
  }

  return true;
}

template<int SKBT>
bool test_set_gradient_at_dual_node();

template<>
bool test_set_gradient_at_dual_node<1>() {
  LCC lcc = create_grid_mesh_with_volume_fraction<6>();
  set_centroids(lcc);

  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point centroid = lcc.attribute<3>(volume)->info().centroid;
    std::array<double, 3> cent = {centroid.x(), centroid.y(), centroid.z()};
    std::array<int, 3> ch = {-1, -1, -1}, flip = {0, 0, 0};
    for(int i = 0; i < 3; i++) {
      if(cent[i] < 1.) {
        ch[i] = 3;
      }
      else if(cent[i] < 2.) {
        ch[i] = 2;
      }
      else if(cent[i] < 3.) {
        ch[i] = 1;
      }
      else if(cent[i] > 10.) {
        ch[i] = 3;
        flip[i] = 1;
      }
      else if(cent[i] > 9.) {
        ch[i] = 2;
        flip[i] = 1;
      }
      else if(cent[i] > 8.) {
        ch[i] = 1;
        flip[i] = 1;
      }
      else if(3. < cent[i] and cent[i] < 8.) {
        ch[i] = 0;
      }
    }
    
    std::array<double, 3> answer;
    if(ch == std::array<int, 3>({0, 0, 0})) {
      answer = {0., 0., 0.};
    }
    else if(ch == std::array<int, 3>({0, 0, 1})) {
      answer = {0., 0., 0.15};
    }
    else if(ch == std::array<int, 3>({0, 1, 0})) {
      answer = {0., 0.15, 0.};
    }
    else if(ch == std::array<int, 3>({1, 0, 0})) {
      answer = {0.15, 0., 0.};
    }
    else if(ch == std::array<int, 3>({0, 0, 2})) {
      answer = {0., 0., 0.5};
    }
    else if(ch == std::array<int, 3>({0, 2, 0})) {
      answer = {0., 0.5, 0.};
    }
    else if(ch == std::array<int, 3>({2, 0, 0})) {
      answer = {0.5, 0., 0.};
    }
    else if(ch == std::array<int, 3>({0, 1, 1})) {
      answer = {0., 0.135, 0.135};
    }
    else if(ch == std::array<int, 3>({1, 0, 1})) {
      answer = {0.135, 0., 0.135};
    }
    else if(ch == std::array<int, 3>({1, 1, 0})) {
      answer = {0.135, 0.135, 0.};
    }
    else if(ch == std::array<int, 3>({0, 1, 2})) {
      answer = {0., 0.085, 0.45};
    }
    else if(ch == std::array<int, 3>({0, 2, 1})) {
      answer = {0., 0.45, 0.085};
    }
    else if(ch == std::array<int, 3>({1, 0, 2})) {
      answer = {0.085, 0., 0.45};
    }
    else if(ch == std::array<int, 3>({1, 2, 0})) {
      answer = {0.085, 0.45, 0.};
    }
    else if(ch == std::array<int, 3>({2, 0, 1})) {
      answer = {0.45, 0., 0.085};
    }
    else if(ch == std::array<int, 3>({2, 1, 0})) {
      answer = {0.45, 0.085, 0.};
    }
    else if(ch == std::array<int, 3>({0, 2, 2})) {
      answer = {0., 1.7/6., 1.7/6.};
    }
    else if(ch == std::array<int, 3>({2, 0, 2})) {
      answer = {1.7/6., 0., 1.7/6.};
    }
    else if(ch == std::array<int, 3>({2, 2, 0})) {
      answer = {1.7/6., 1.7/6., 0.};
    }
    else if(ch == std::array<int, 3>({1, 1, 1})) {
      answer = {2.187/18., 2.187/18., 2.187/18.}; // 計算間違いかも
    }
    else if(ch == std::array<int, 3>({1, 1, 2})) {
      answer = {0.0765, 0.0765, 0.405};
    }
    else if(ch == std::array<int, 3>({1, 2, 1})) {
      answer = {0.0765, 0.405, 0.0765};
    }
    else if(ch == std::array<int, 3>({2, 1, 1})) {
      answer = {0.405, 0.0765, 0.0765};
    }
    else if(ch == std::array<int, 3>({1, 2, 2})) {
      answer = {0.867/18., 0.255, 0.255};
    }
    else if(ch == std::array<int, 3>({2, 1, 2})) {
      answer = {0.255, 0.867/18., 0.255};
    }
    else if(ch == std::array<int, 3>({2, 2, 1})) {
      answer = {0.255, 0.255, 0.867/18.};
    }
    else if(ch == std::array<int, 3>({2, 2, 2})) {
      answer = {2.89/18., 2.89/18., 2.89/18.};
    }
    else if(ch[0] == 3 or ch[1] == 3 or ch[2] == 3) {
      continue;
    }
    else {
      assert(false);
    }

    for(int i = 0; i < 3; i++) {
      if(flip[i]) answer[i] = -answer[i];
    }

    Vector answer_vec = {answer[0], answer[1], answer[2]};

    __set_gradient_at_dual_node(lcc, volume);
    Vector result_vec = lcc.attribute<3>(volume)->info().gradient;
    if(!is_approximately_equal_Vector(answer_vec, result_vec)) {
      // std::cout << "CENTER: " << lcc.attribute<3>(volume)->info().centroid << "\nANSWER: " << answer_vec << "\nRESULT: " << result_vec << '\n' << std::endl;
      return false;
    }
  }
  return true;
}

template<int SKBT>
std::vector<Point> test_laplacian_smoothing_for_unmarked_cells();
template<int SKBT>
std::vector<Point> answer_laplacian_smoothing_for_unmarked_cells();

template<>
std::vector<Point> test_laplacian_smoothing_for_unmarked_cells<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  size_type no_mark = lcc.get_new_mark();
  auto vertices = lcc.one_dart_per_cell<0>();

  int cnt = 0;
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    if(cnt++&1) {
      lcc.mark_cell<0>(vertex, no_mark);
    }
  }

  laplacian_smoothing_for_unmarked_cells(lcc, no_mark);

  std::vector<Point> result;
  result.reserve(125);

  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    result.emplace_back(lcc.point(vertex));
  }

  return result;
}

template<>
std::vector<Point> answer_laplacian_smoothing_for_unmarked_cells<1>() {
  std::vector<Point> answer(125);
  for(int i = 0; i < 5; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) {
    answer[i*25+j*5+k] = {double(i), double(j), double(k)};
  }
  return answer;
}

template<>
std::vector<Point> test_laplacian_smoothing_for_unmarked_cells<2>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  size_type l_mark = lcc.get_new_mark();
  auto vertices = lcc.one_dart_per_cell<0>();

  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    int x = p.x()+0.1, y = p.y()+0.1, z = p.z()+0.1;
    if(x == 2 and y == 2) {
      if(z == 1) {
        lcc.point(vertex) = {2, 2, 3};
      }
      else if(z == 3) {
        lcc.point(vertex) = {2, 2, 1};
      }
    }
    else if(x == 2 and z == 2) {
      if(y == 1) {
        lcc.point(vertex) = {2, 3, 2};
      }
      else if(y == 3) {
        lcc.point(vertex) = {2, 1, 2};
      }
    }
    else if(y == 2 and z == 2) {
      if(x == 1) {
        lcc.point(vertex) = {3, 2, 2};
      }
      else if(x == 3) {
        lcc.point(vertex) = {1, 2, 2};
      }
    }
    if(z > 2) lcc.mark_cell<0>(vertex, l_mark);
  }

  laplacian_smoothing_for_unmarked_cells(lcc, l_mark);

  std::vector<Point> result;
  result.reserve(125);

  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    result.emplace_back(lcc.point(vertex));
    // std::cout << lcc.point(vertex) << std::endl;
  }

  return result;
}

template<>
std::vector<Point> answer_laplacian_smoothing_for_unmarked_cells<2>() {
  std::vector<Point> answer = {{4./3, 4./3, 2}, {4./3, 2, 4./3}, {2, 4./3, 4./3}, {8./3, 4./3, 2}, {8./3, 2, 4./3}, {2, 8./3, 4./3}, {4./3, 8./3, 2}, {8./3, 8./3, 2}, {2, 2, 1}};
  answer.reserve(125);

  for(int i = 0; i < 5; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) {
    if(i == 0 or i == 4 or j == 0 or j == 4 or k == 0 or k == 4) {
      answer.emplace_back(i, j, k);
    }
    else if(((i^j^k)&1) or (i == 2 and j == 2 and k == 2) or k > 2) {
      if(i == 2 and j == 2 and k == 3) continue;
      answer.emplace_back(i, j, k);
    }
  }

  return answer;
}

template<int SKBT>
bool test2_laplacian_smoothing_for_unmarked_cells();

template<>
bool test2_laplacian_smoothing_for_unmarked_cells<1>() {
  LCC lcc = create_grid_mesh(5, 5, 5);
  size_type no_mark = lcc.get_new_mark();
  auto vertices = lcc.one_dart_per_cell<0>();
  set_centroids(lcc);

  const Point C = {2.5, 2.5, 2.5};
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point& p = lcc.point(vertex);
    for(int i = 1; i < 5; i++) for(int j = 1; j < 5; j++) for(int k = 1; k < 5; k++) {
      if(!is_approximately_equal_Point(p, {i, j, k})) continue;
      if(i == 1 or i == 4 or j == 1 or j == 4 or k == 1 or k == 4) {
        p = C + 0.3*(p-C);
      }
      else {
        p = C + 7.*(p-C);
      }
    }
  }

  int inv_counter = 0;
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    std::array<typename LCC::Vertex_attribute_descriptor, 8> nodes = {lcc.attribute<0>(volume), lcc.attribute<0>(lcc.beta(volume, 1)), lcc.attribute<0>(lcc.beta(volume, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 0)), 
                                                                      lcc.attribute<0>(lcc.beta(volume, 0, 2, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 2, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 2, 0)), lcc.attribute<0>(lcc.beta(volume, 1, 2, 0))};
    double j = scale_jacobian_of_hexa(lcc, nodes);
    if(j < 0.) inv_counter++;
  }
  if(inv_counter != 26) return false;

  laplacian_smoothing_for_unmarked_cells(lcc, no_mark);
  laplacian_smoothing_for_unmarked_cells(lcc, no_mark);

  inv_counter = 0;
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    std::array<typename LCC::Vertex_attribute_descriptor, 8> nodes = {lcc.attribute<0>(volume), lcc.attribute<0>(lcc.beta(volume, 1)), lcc.attribute<0>(lcc.beta(volume, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 0)), 
                                                                      lcc.attribute<0>(lcc.beta(volume, 0, 2, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 2, 1, 1)), lcc.attribute<0>(lcc.beta(volume, 2, 0)), lcc.attribute<0>(lcc.beta(volume, 1, 2, 0))};
    double j = scale_jacobian_of_hexa(lcc, nodes);
    if(j < 0.) inv_counter++;
  }

  return inv_counter == 0;
}

template<int SKBT>
bool test_move_points_onto_mesh();

// moving points onto plane x+y+z=5.1
template<>
bool test_move_points_onto_mesh<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  const double xyz_sum = 5.1;
  size_type move_mark = lcc.get_new_mark();
  auto detectIntersection = [&](LCC& lcc, Dart_handle dart) -> bool {
    auto &face_attr = lcc.attribute<2>(dart)->info();
    Point p = face_attr.dual_edge.source();
    Point q = face_attr.dual_edge.target();
    double p_sum = p.x()+p.y()+p.z();
    double q_sum = q.x()+q.y()+q.z();
    if((p_sum < xyz_sum and q_sum < xyz_sum) or (p_sum > xyz_sum and q_sum > xyz_sum)) {
      return false;
    }

    Vector p_vec = p-CGAL::ORIGIN;
    Vector q_vec = q-CGAL::ORIGIN;
    
    face_attr.intersection = CGAL::ORIGIN + p_vec + ((xyz_sum-p_sum)/(q_sum-p_sum))*(q_vec-p_vec);
    face_attr.normal = {1./sqrt(3.), 1./sqrt(3.), 1./sqrt(3.)};

    return true;
  };
  move_points_onto_mesh(lcc, move_mark, detectIntersection);

  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    
    double x = p.x(), y = p.y(), z = p.z();
    if(!is_approximately_equal_Point(p, Point(int(x+0.1), int(y+0.1), int(z+0.1)))) {
      if(!lcc.is_marked(vertex, move_mark)) {
        return false;
      }

      if(!is_approximately_equal_double(x+y+z, xyz_sum)) {
        return false;
      }
      double dx = x-int(x+100)-100, dy = y-int(y+100)-100, dz = z-int(z+100)-100;
      if(!is_approximately_equal_double(dx, dy) or !is_approximately_equal_double(dy, dz) or !is_approximately_equal_double(dx, dz)) {
        return false;
      }
    }
    else {
      if(lcc.is_marked(vertex, move_mark)) {
        return false;
      }
    }
  }
  
  return true;
}

// moving points onto plane x=1.3 (y < 1.3) and y=1.3 (x < 1.3)
template<>
bool test_move_points_onto_mesh<2>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  const double xy_size = 1.3;
  size_type move_mark = lcc.get_new_mark();
  auto detectIntersection = [&](LCC& lcc, Dart_handle dart) -> bool {
    auto &face_attr = lcc.attribute<2>(dart)->info();
    Point p = face_attr.dual_edge.source();
    Point q = face_attr.dual_edge.target();
    if((p.x() > xy_size or p.y() > xy_size) and (q.x() > xy_size or q.y() > xy_size)) {
      return false;
    }
    if((p.x() < xy_size and p.y() < xy_size) and (q.x() < xy_size and q.y() < xy_size)) {
      return false;
    }

    Vector p_vec = p-CGAL::ORIGIN;
    Vector q_vec = q-CGAL::ORIGIN;

    if(abs(p.x()-q.x()) > 0.5) {
      face_attr.intersection = {xy_size, p.y(), p.z()};
      face_attr.normal = {1., 0., 0.};
    }
    else if(abs(p.y()-q.y()) > 0.5) {
      face_attr.intersection = {p.x(), xy_size, p.z()};
      face_attr.normal = {0., 1., 0.};
    }

    return true;
  };
  move_points_onto_mesh(lcc, move_mark, detectIntersection);

  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    
    double x = p.x(), y = p.y(), z = p.z();
    if(!is_approximately_equal_Point(p, Point(int(x+0.1), int(y+0.1), int(z+0.1)))) {
      if(!lcc.is_marked(vertex, move_mark)) {
        return false;
      }
      
      if(!is_approximately_equal_Point(p, {(1.+xy_size)*0.5, (1.+xy_size)*0.5, p.z()}) and !is_approximately_equal_Point(p, {0., xy_size, p.z()}) and !is_approximately_equal_Point(p, {xy_size, 0., p.z()})) {
        return false;
      }
    }
    else {
      if(lcc.is_marked(vertex, move_mark)) {
        return false;
      }
    }
  }
  
  return true;
}

template<int>
bool test_volumes_around_node();

template<>
bool test_volumes_around_node<1>() {
  LCC lcc = create_grid_mesh(2, 2, 2);
  set_centroids(lcc);
  std::array<int, 6> available_codes = {0b000'100'110'010'001'101'111'011,
                                        0b000'010'110'100'001'011'111'101,
                                        0b000'001'101'100'010'011'111'110,
                                        0b000'100'101'001'010'110'111'011,
                                        0b000'010'011'001'100'110'111'101,
                                        0b000'001'011'010'100'101'111'110};
  auto checker = [&](std::array<Dart_handle, 8> volumes)->bool {
    Point base = lcc.attribute<3>(volumes[0])->info().centroid;
    int code = 0, base_code = !is_approximately_equal_double(base.x(), 0.5) ^ !is_approximately_equal_double(base.y(), 0.5) ^ !is_approximately_equal_double(base.z(), 0.5);
    for(int i = 1; i < 8; i++) {
      Point now = lcc.attribute<3>(volumes[i])->info().centroid;
      code <<= 3;
      code += !is_approximately_equal_double(base.x(), now.x())*4 + !is_approximately_equal_double(base.y(), now.y())*2 + !is_approximately_equal_double(base.z(), now.z());
    }
    for(int i = 0; i < 6; i++) {
      if(available_codes[i] != code) continue;
      if(base_code^(i&1)^1) {
        return true;
      }
      // std::cout << "id: " << i << std::endl;
      // std::cout << "FAILURE: " << lcc.attribute<3>(volumes[0])->info().centroid << ' ' << lcc.point(lcc.beta<1>(volumes[0])) << ' ' << lcc.point(lcc.beta(volumes[0], 1, 1)) << std::endl;
    }
    return false;
  };
  
  auto darts = lcc.darts();
  for(auto dart = darts.begin(); dart != darts.end(); dart++) {
    if(is_approximately_equal_Point({1, 1, 1}, lcc.point(dart))) {
      std::array<Dart_handle, 8> volumes = volumes_around_node(lcc, dart);
      if(!checker(volumes)) {
        return false;
      }
    }
  }
  return true;
}

template<int>
bool test_get_signal();

template<>
bool test_get_signal<1>() {
  LCC lcc = create_grid_mesh(2, 2, 2);
  size_type test_mark = lcc.get_new_mark();
  std::array<Dart_handle, 8> volumes;
  set_centroids(lcc); // just for setting attributes on 3-cells

  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    auto now_volumes = volumes_around_node(lcc, vertex);
    bool is_invalid = false;
    for(auto volume: now_volumes) {
      if(volume == nullptr or lcc.attribute<3>(volume) == nullptr)
        is_invalid = true;
    }
    
    if(is_invalid) {
      if(__get_signal(lcc, vertex, test_mark) != -1) {
        return false;
      }
    }
    else {
      volumes = now_volumes;
    }
  }

  for(int i = 0; i < (1<<8); i++) {
    for(int j = 0; j < 8; j++) if(i>>j&1) {
      lcc.mark_cell<3>(volumes[7-j], test_mark);
    }

    if(__get_signal(lcc, volumes[0], test_mark) != i) {
      return false;
    }
    
    lcc.unmark_all(test_mark);
  }
  lcc.free_mark(test_mark);

  return true;
}

template<int>
bool test_get_seven_non_manifold_templates();

template<>
bool test_get_seven_non_manifold_templates<1>() {
  std::array<std::vector<int>, 7> seven_templates = __get_seven_non_manifold_templates();
  std::array<int, 256> templates = __get_non_manifold_template_list();
  for(int i = 0; i < 256; i++) {
    bool flag = false;
    for(int j = 7; j--;) for(auto temp: seven_templates[j]) {
      if(i == temp) flag = true;
    }
    if(templates[i]) {
      if(!flag) return false;
    }
    else {
      if(flag) return false;
    }
  }

  return true;
}

template<int>
bool test_get_non_manifold_template_list();

template<>
bool test_get_non_manifold_template_list<1>() {
  std::array<int, 256> templates = __get_non_manifold_template_list();
  for(int i = 0; i < 256; i++) {
    std::array<std::array<bool, 8>, 8> graph = {};
    graph[0][1] = graph[1][0] = !((i^(i>>1))&1);
    graph[0][3] = graph[3][0] = !((i^(i>>3))&1);
    graph[0][4] = graph[4][0] = !((i^(i>>4))&1);
    graph[1][2] = graph[2][1] = !((i^(i>>1))&2);
    graph[1][5] = graph[5][1] = !((i^(i>>4))&2);
    graph[2][3] = graph[3][2] = !((i^(i>>1))&4);
    graph[2][6] = graph[6][2] = !((i^(i>>4))&4);
    graph[3][7] = graph[7][3] = !((i^(i>>4))&8);
    graph[4][5] = graph[5][4] = !((i^(i>>1))&16);
    graph[4][7] = graph[7][4] = !((i^(i>>3))&16);
    graph[5][6] = graph[6][5] = !((i^(i>>1))&32);
    graph[6][7] = graph[7][6] = !((i^(i>>1))&64);
    
    std::array<int, 8> group = {-1, -1, -1, -1, -1, -1, -1, -1};
    int group_id = 0;
    for(int j = 0; j < 8; j++) if(group[j] == -1) {
      std::queue<int> q;
      q.emplace(j);
      while(!q.empty()) {
        int now = q.front();
        q.pop();

        if(group[now] != -1) continue;
        group[now] = group_id;
        for(int k = 0; k < 8; k++) {
          if(graph[now][k]) q.emplace(k);
        }
      }

      group_id++;
    }
    if(int(group_id > 2) ^ templates[i]) {
      return false;
    }
  }

  return true;
}

template<int>
bool test_get_solution_to_non_manifold_templates_list();

template<>
bool test_get_solution_to_non_manifold_templates_list<1>() {
  std::array<int, 256> is_in_templates = __get_non_manifold_template_list();
  std::array<std::vector<int>, 256> resolve_templates = __get_solution_to_non_manifold_templates_list();
  
  for(int i = 0; i < 256; i++) {
    if(!is_in_templates[i]) {
      if(!resolve_templates[i].empty()) {
        return false;
      }
      continue;
    }

    // check if the solutions are ordered in bit-popcount
    int pp_counter = 0;
    for(int sol: resolve_templates[i]) {
      int now = __builtin_popcount(sol);
      if(pp_counter > now) {
        return false;
      }
      pp_counter = now;
    }

    // check if the solutions are correct and complete
    std::array<int, 256> sol_list = {};
    for(int sol: resolve_templates[i]) {
      sol_list[sol] = 1;
    }
    for(int j = 0; j < 256; j++) {
      int now = i^j;
      if(sol_list[j] ^ !is_in_templates[now]) {
        return false;
      }
    }
  }

  return true;
}

template<int>
bool test_resolve_non_manifold_case();

// x+y+z < 5.1, thus no non-manifold case
template<>
bool test_resolve_non_manifold_case<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  auto volumes = lcc.one_dart_per_cell<3>();
  set_centroids(lcc);

  const double c_xyz = 5.1;

  size_type inner_mark = lcc.get_new_mark();
  std::vector<double> fracs;

  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point p = lcc.attribute<3>(volume)->info().centroid;
    double &frac = lcc.attribute<3>(volume)->info().fraction;
    double x = p.x(), y = p.y(), z = p.z();
    double sum_xyz = x+y+z;
    if(sum_xyz < c_xyz-1.5) {
      frac = 1.;
    }
    else if(sum_xyz < c_xyz-0.5) {
      double l = int(c_xyz)+1. - c_xyz;
      frac = 1. - l*l*l/6.;
    }
    else if(sum_xyz < c_xyz+0.5) {
      double l = c_xyz - int(c_xyz);
      frac = (l+1)*(l+1)*(l+1)/6. - l*l*l/2.;
    }
    else if(sum_xyz < c_xyz+1.5) {
      double l = c_xyz - int(c_xyz);
      frac = l*l*l/6.;
    }
    else {
      frac = 0.;
    }

    if(sum_xyz < c_xyz) {
      lcc.mark_cell<3>(volume, inner_mark);
    }
    fracs.emplace_back(frac);
  }

  resolve_non_manifold_case(lcc, 0.5, inner_mark);

  std::vector<double> after_fracs;
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    double frac = lcc.attribute<3>(volume)->info().fraction;
    after_fracs.emplace_back(lcc.attribute<3>(volume)->info().fraction);
    if((frac <= 0.5 and lcc.is_marked(volume, inner_mark)) or (frac > 0.5 and !lcc.is_marked(volume, inner_mark))) {
      return false;
    }
  }

  for(int i = 0; i < 64; i++) {
    if(!is_approximately_equal_double(after_fracs[i], fracs[i])) {
      return false;
    }
  }

  lcc.free_mark(inner_mark);

  return true;
}

// x+y+z < 4.8 and x+y+z > 7.2
template<>
bool test_resolve_non_manifold_case<2>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  auto volumes = lcc.one_dart_per_cell<3>();
  set_centroids(lcc);

  const double c_xyz = 4.8;
  const double d_xyz = 7.2;

  size_type inner_mark = lcc.get_new_mark();

  set_plane(lcc, c_xyz, d_xyz, inner_mark);
  
  resolve_non_manifold_case(lcc, 0.5, inner_mark);
  
  auto is_in_template = __get_non_manifold_template_list();
  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    int signal = __get_signal(lcc, vertex, inner_mark);
    if(signal == -1) continue;
    if(is_in_template[signal]) {
      return false;
    }
  }

  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
    if(vol_attr.fraction < 0.-1e-10 or vol_attr.fraction > 1.+1e-10) {
      return false;
    }
    if((vol_attr.fraction <= 0.5 and lcc.is_marked(volume, inner_mark)) or (vol_attr.fraction > 0.5 and !lcc.is_marked(volume, inner_mark))) {
      return false;
    }
  }

  lcc.free_mark(inner_mark);

  return true;
}

// x+y+z < 5.4 and x+y+z > 6.6
template<>
bool test_resolve_non_manifold_case<3>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  auto volumes = lcc.one_dart_per_cell<3>();
  set_centroids(lcc);

  const double c_xyz = 5.4;
  const double d_xyz = 6.6;

  size_type inner_mark = lcc.get_new_mark();

  set_plane(lcc, c_xyz, d_xyz, inner_mark);
  
  resolve_non_manifold_case(lcc, 0.5, inner_mark);
  
  auto is_in_template = __get_non_manifold_template_list();
  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    int signal = __get_signal(lcc, vertex, inner_mark);
    if(signal == -1) continue;
    if(is_in_template[signal]) {
      return false;
    }
  }

  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
    if(vol_attr.fraction < 0.-1e-10 or vol_attr.fraction > 1.+1e-10) {
      return false;
    }
    if((vol_attr.fraction <= 0.5 and lcc.is_marked(volume, inner_mark)) or (vol_attr.fraction > 0.5 and !lcc.is_marked(volume, inner_mark))) {
      return false;
    }
  }

  lcc.free_mark(inner_mark);

  return true;
}

// x+y+z < 5.4 and x+y+z > 6.6 and a little bit change
// this will fail because it violate the prerequisite of resolve_non_manifold_case
template<>
bool test_resolve_non_manifold_case<4>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  auto volumes = lcc.one_dart_per_cell<3>();
  set_centroids(lcc);

  const double c_xyz = 5.4;
  const double d_xyz = 6.6;

  size_type inner_mark = lcc.get_new_mark();

  set_plane(lcc, c_xyz, d_xyz, inner_mark);
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
    if(vol_attr.centroid.x()+vol_attr.centroid.y()+vol_attr.centroid.z() > 6. and vol_attr.fraction < 0.5-1e-10) {
      vol_attr.fraction = 0.1;
    }
  }
  
  resolve_non_manifold_case(lcc, 0.5, inner_mark);
  
  auto is_in_template = __get_non_manifold_template_list();
  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    int signal = __get_signal(lcc, vertex, inner_mark);
    if(signal == -1) continue;
    if(is_in_template[signal]) {
      return false;
    }
  }

  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
    if(vol_attr.fraction < 0.-1e-10 or vol_attr.fraction > 1.+1e-10) {
      return false;
    }
    if((vol_attr.fraction <= 0.5 and lcc.is_marked(volume, inner_mark)) or (vol_attr.fraction > 0.5 and !lcc.is_marked(volume, inner_mark))) {
      return false;
    }
  }

  lcc.free_mark(inner_mark);

  return true;
}

template<int>
bool test_detect_intersection_with_volume_fraction();

template<>
bool test_detect_intersection_with_volume_fraction<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  const double xyz_sum = 5.1;
  size_type move_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();

  set_plane(lcc, xyz_sum, 100., inner_mark);
  
  auto detectIntersection = detect_intersection_with_volume_fraction(0.5, inner_mark, move_mark);

  set_dual_edges(lcc);

  auto faces = lcc.one_dart_per_cell<2>();
  for(auto face = faces.begin(); face != faces.end(); face++) {
    auto &face_attr = lcc.attribute<2>(face)->info();
    Segment dual_edge = face_attr.dual_edge;
    Point p = dual_edge.source(), q = dual_edge.target();
    double px = p.x(), py = p.y(), pz = p.z();
    double qx = q.x(), qy = q.y(), qz = q.z();
    
    if(detectIntersection(lcc, face)) {
      if(lcc.is_free<3>(face)) {
        return false;
      }
      if(!(int(lcc.is_marked(face, inner_mark)) ^ int(lcc.is_marked(lcc.beta<3>(face), inner_mark)))) {
        return false;
      }
      double frac1 = lcc.attribute<3>(face)->info().fraction, frac2 = lcc.attribute<3>(lcc.beta<3>(face))->info().fraction;
      Point inter = face_attr.intersection;
      Vector norm = face_attr.normal;
      if(!is_approximately_equal_Point(inter, p + ((0.5-frac1)/(frac2-frac1))*(q-p))) {
        return false;
      }
      if(1<px and 1<py and 1<pz and px<3 and py<3 and pz<3 and 1<qx and 1<qy and 1<qz and qx<3 and qy<3 and qz<3) {
        if(!is_approximately_equal_Vector(norm, {-1./sqrt(3.), -1./sqrt(3.), -1./sqrt(3.)})) {
          return false;
        }
      }
    //   std::cout << "fraction 1: " << lcc.attribute<3>(face)->info().fraction << " : " << p << std::endl;
    //   std::cout << "fraction 2: " << lcc.attribute<3>(lcc.beta<3>(face))->info().fraction << " : " << q << std::endl;
    //   std::cout << "intersection: " << inter << std::endl;
    //   std::cout << "normal: " << norm << std::endl;
    }
    else {
      if(!lcc.is_free<3>(face) and (int(lcc.is_marked(face, inner_mark)) ^ int(lcc.is_marked(lcc.beta<3>(face), inner_mark)))) {
        return false;
      }
    }
  }

  return true;
}

template<int>
bool test_move_points_onto_mesh_with_volume_fraction();

// moving points onto plane x+y+z=5.1
// (2, 2, 2) should go to (a, a, a) where a is similar to 1.7.
template<>
bool test_move_points_onto_mesh_with_volume_fraction<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  const double xyz_sum = 5.1;
  size_type no_mark = lcc.get_new_mark();
  size_type move_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();

  set_plane(lcc, xyz_sum, 100., no_mark);

  auto volumes = lcc.one_dart_per_cell<3>();
  double frac1, frac2;
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
    if(is_approximately_equal_Point(vol_attr.centroid, {1.5, 2.5, 0.5})) {
      frac1 = vol_attr.fraction;
    }
    if(is_approximately_equal_Point(vol_attr.centroid, {3.5, 1.5, 0.5})) {
      frac2 = vol_attr.fraction;
    }
  }
  
  move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    auto &vol_attr = lcc.attribute<3>(volume)->info();
  }

  auto vertices = lcc.one_dart_per_cell<0>();
  int count = 0;
  double coord_sum = int(xyz_sum) - 0.5 + (0.5-frac1)/(frac2-frac1);
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    if(is_approximately_equal_Point(p, {coord_sum/3., coord_sum/3., coord_sum/3.})) {
      count++;
    }
    if(int(!is_approximately_equal_Point(p, {int(p.x()+1e-10), int(p.y()+1e-10), int(p.z()+1e-10)})) ^ int(lcc.is_marked(vertex, move_mark))) {
      return false;
    }
  }
  if(count != 1) {
    return false;
  }

  lcc.free_mark(move_mark);
  lcc.free_mark(inner_mark);

  return true;
}

// Since move_points_onto_mesh_with_volume_fraction approximates mesh and doesn't use mesh itself, it is difficult to estimate the exact destination of each node.


template<int>
bool test_get_orthogonal_vectors();

template<>
bool test_get_orthogonal_vectors<1>() {
  auto check = [](Vector original)->bool {
    auto [out1, out2] = get_orthogonal_vectors(original);
    if(!is_approximately_equal_double(out1*out1, 1.) or !is_approximately_equal_double(out2*out2, 1.)) {
      return false;
    }
    if(!is_approximately_equal_double(original*out1, 0.) or !is_approximately_equal_double(original*out2, 0.)) {
      return false;
    }
    if(!is_approximately_equal_double(out1*out2, 0.)) {
      return false;
    }
    return true;
  };

  if(!check({1., 0., 0.})) {
    return false;
  }
  if(!check({0., 1., 0.})) {
    return false;
  }
  if(!check({0., 0., 1.})) {
    return false;
  }

  RandomPointGenerator gen(2.);
  int number_of_random_points = 100;
  std::vector<Point> random_points;
  copy_n(gen, number_of_random_points, std::back_inserter(random_points));

  for(auto point: random_points) {
    if(!check(point-CGAL::ORIGIN)) {
      return false;
    }
  }

  return true;
}

template<int>
bool test_get_neighbors_list_for_smoothing();

template<>
bool test_get_neighbors_list_for_smoothing<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  size_type no_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();
  size_type move_mark = lcc.get_new_mark();
  set_plane(lcc, 5.1, 100, no_mark);
  // CGAL::HexRefinement::render_two_refinement_result_with_mark(lcc, no_mark);
  move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);

  std::vector<std::vector<Dart_handle>> neighbors_list = get_neighbors_list_for_smoothing(lcc, move_mark, inner_mark);

  lcc.free_mark(no_mark);
  lcc.free_mark(inner_mark);
  lcc.free_mark(move_mark);

  int cnt = 0;
  std::vector<int> pos(125, -1);
  for(int i = 0; i < 125; i++) {
    if(neighbors_list[i].empty()) continue;
    pos[i] = cnt++;
  }

  std::vector<std::pair<int, int>> edges_result;
  for(int i = 0; i < 125; i++) {
    for(auto neighbor: neighbors_list[i]) {
      int id = lcc.attribute<0>(neighbor)->id;
      if(i < id) continue;
      edges_result.emplace_back(pos[i], pos[id]);
    }
  }
  std::vector<std::pair<int, int>> edges_ideal = {
    {0, 1}, {1, 2}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 8}, {8, 5},
    {5, 9}, {9, 7}, {3, 11}, {8, 12}, {5, 13}, {9, 14}, {7, 15}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
    {14, 15}, {15, 16}, {10, 21}, {10, 17}, {17, 22}, {17, 12}, {12, 23}, {12, 18}, {18, 24}, {18, 14}, {14, 25},
    {14, 19}, {19, 26}, {19, 16}, {16, 27}, {20, 21}, {21, 22}, {22, 23}, {23, 24}, {24, 25}, {25, 26}, {26, 27},
    {27, 28}, {20, 29}, {29, 33}, {29, 22}, {22, 34}, {22, 30}, {30, 35}, {30, 24}, {24, 36}, {24, 31}, {31, 37},
    {31, 26}, {26, 38}, {26, 32}, {32, 39}, {32, 28}, {33, 34}, {34, 35}, {35, 36}, {36, 37}, {37, 38}, {38, 39}
  };
  const boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> graph_result(edges_result.begin(), edges_result.end(), cnt), graph_ideal(edges_ideal.begin(), edges_ideal.end(), 40);

  return boost::isomorphism(graph_result, graph_ideal);
}

template<>
bool test_get_neighbors_list_for_smoothing<2>() {
  LCC lcc = create_grid_mesh(4, 4, 4);
  set_centroids(lcc);
  size_type no_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();
  size_type move_mark = lcc.get_new_mark();

  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point p = lcc.attribute<3>(volume)->info().centroid;
    double x = p.x(), y = p.y(), z = p.z();
    if(x < 1. or y < 1. or (x < 3. and y < 3. and 1. < z and z < 2.)) {
      lcc.attribute<3>(volume)->info().fraction = 0.8;
      lcc.mark_cell<3>(volume, no_mark);
    }
    else {
      lcc.attribute<3>(volume)->info().fraction = 0.;
    }
  }
  // CGAL::HexRefinement::render_two_refinement_result_with_mark(lcc, no_mark);
  move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);

  std::vector<std::vector<Dart_handle>> neighbors_list = get_neighbors_list_for_smoothing(lcc, move_mark, inner_mark);

  lcc.free_mark(no_mark);
  lcc.free_mark(inner_mark);
  lcc.free_mark(move_mark);

  int cnt = 0;
  std::vector<int> pos(125, -1);
  for(int i = 0; i < 125; i++) {
    if(neighbors_list[i].empty()) continue;
    pos[i] = cnt++;
  }

  std::vector<std::pair<int, int>> edges_result;
  for(int i = 0; i < 125; i++) {
    for(auto neighbor: neighbors_list[i]) {
      int id = lcc.attribute<0>(neighbor)->id;
      if(i < id) continue;
      edges_result.emplace_back(pos[i], pos[id]);
    }
  }
  std::vector<std::pair<int, int>> edges_ideal = {
    {0, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {5, 10}, {6, 11}, {7, 12}, {8, 13}, {9, 14},
    {10, 11}, {11, 12}, {13, 14}, {10, 15}, {11, 16}, {12, 17}, {13, 18}, {14, 19}, {15, 16}, {16, 17}, {18, 19},
    {15, 20}, {16, 21}, {17, 22}, {18, 23}, {19, 24}, {20, 21}, {21, 22}, {23, 24}, {20, 25}, {21, 26}, {22, 27},
    {23, 28}, {24, 29}, {25, 26}, {26, 27}, {27, 28}, {28, 29}, {25, 30}, {26, 31}, {27, 32}, {28, 33}, {29, 34},
    {30, 31}, {31, 32}, {32, 33}, {33, 34}, {7, 39}, {8, 40}, {12, 35}, {13, 36}, {39, 41}, {40, 42}, {35, 37}, {36, 38},
    {22, 35}, {35, 39}, {27, 37}, {37, 41}, {23, 36}, {36, 40}, {28, 38}, {38, 42}, {39, 40}, {41, 42}, {37, 38}
  };
  const boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> graph_result(edges_result.begin(), edges_result.end(), cnt), graph_ideal(edges_ideal.begin(), edges_ideal.end(), 40);

  return boost::isomorphism(graph_result, graph_ideal);
}

template<int>
bool test_surface_smoothing();

// see laplacian smoothing on a plane 0.1x+0.1y+0.8z = 2
template<>
bool test_surface_smoothing<1>() {
  LCC lcc = create_grid_mesh(4, 4, 4);

  std::array<std::array<Dart_handle, 5>, 5> darts;

  size_type surface_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();

  set_centroids(lcc);
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point p = lcc.attribute<3>(volume)->info().centroid;
    if(p.z() < 2.) {
      lcc.mark_cell<3>(volume, inner_mark);
    }
  }

  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    double x = p.x(), y = p.y(), z = p.z();
    if(z < 1.5 or 2.5 < z) continue;
    for(int i = 0; i < 5; i++) for(int j = 0; j < 5; j++) {
      if(is_approximately_equal_Point(p, {i, j, 2.})) {
        darts[i][j] = vertex;
      }
    }
    if(is_approximately_equal_Point(p, {2., 2., 2.})) {
      x = 1.5;
      y = 1.5;
    }
    lcc.point(vertex) = {x, y, 2.5 - 0.125*(x+y)};

    lcc.attribute<0>(vertex)->normal = {1./sqrt(66.), 1./sqrt(66.), 8./sqrt(66.)};

    lcc.mark_cell<0>(vertex, surface_mark);
  }

  surface_smoothing(lcc, surface_mark, inner_mark);

  lcc.free_mark(surface_mark);
  lcc.free_mark(inner_mark);

  std::array<std::array<std::pair<double, double>, 5>, 5> answer = {{{{{.5, .5}, {1./3, 1.}, {1./3, 2.}, {1./3, 3.}, {.5, 3.5}}},
                                                                      {{{1., 1./3}, {1., 1.}, {3.5/4, 7.5/4}, {1., 3.}, {1., 11./3}}},
                                                                      {{{2., 1./3}, {7.5/4, 3.5/4}, {2., 2.}, {7.5/4, 11.5/4}, {2., 11./3}}},
                                                                      {{{3., 1./3}, {3., 1.}, {11.5/4, 7.5/4}, {3., 3.}, {3., 11./3}}},
                                                                      {{{3.5, .5}, {11./3, 1.}, {11./3, 2.}, {11./3, 3.}, {3.5, 3.5}}}}};

  for(int i = 0; i < 5; i++) for(int j = 0; j < 5; j++) {
    auto [x, y] = answer[i][j];
    if(!is_approximately_equal_Point(lcc.point(darts[i][j]), {x, y, 2.5 - 0.125*(x+y)})) {
      return false;
    }
  }

  return true;
}

// see laplacian smoothing on a curved surface z = 2 + 0.1(x-2)^2 - 0.1(y-2)^2 with normal (0, 0, 1)
template<>
bool test_surface_smoothing<2>() {
  LCC lcc = create_grid_mesh(4, 4, 4);

  std::array<std::array<Dart_handle, 5>, 5> darts;

  size_type surface_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();

  set_centroids(lcc);
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    Point p = lcc.attribute<3>(volume)->info().centroid;
    if(p.z() < 2.) {
      lcc.mark_cell<3>(volume, inner_mark);
    }
  }

  auto vertices = lcc.one_dart_per_cell<0>();
  for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
    Point p = lcc.point(vertex);
    double x = p.x(), y = p.y(), z = p.z();
    if(z < 1.5 or 2.5 < z) continue;
    for(int i = 0; i < 5; i++) for(int j = 0; j < 5; j++) {
      if(is_approximately_equal_Point(p, {i, j, 2.})) {
        darts[i][j] = vertex;
      }
    }
    if(is_approximately_equal_double(p.x(), 1.)) {
      x = 1.5;
    }
    if(is_approximately_equal_double(p.y(), 1.)) {
      y = 1.5;
    }
    lcc.point(vertex) = {x, y, 0.1*(x*x - y*y) - 0.4*(x - y) + 2.};

    lcc.attribute<0>(vertex)->normal = {0., 0., 1.};

    lcc.mark_cell<0>(vertex, surface_mark);
  }

  surface_smoothing(lcc, surface_mark, inner_mark, 0.);

  lcc.free_mark(surface_mark);
  lcc.free_mark(inner_mark);

  // cannot assure the destination of the points i=0,5, j=0,5, because there are not enough information for the matrices to have rank 4.
  for(int i = 1; i < 4; i++) for(int j = 1; j < 4; j++) {
    // std::cout << "{" << i << ", " << j << "} : " << lcc.point(darts[i][j]) << std::endl;
    Point p = lcc.point(darts[i][j]);
    if(!is_approximately_equal_double(p.z(), 0.1*(p.x()*p.x()-p.y()*p.y()) - 0.4*(p.x() - p.y()) + 2.)) {
      return false;
    }
  }

  return true;
}



int main() {
  
  TestFramework test;

  // テストいっぱい
  test.test("test set_dual_edges", test_set_dual_edges<1>);
  test.test("test __set_gradient_at_dual_node", test_set_gradient_at_dual_node<1>);
  test.test_approximate_value<std::vector<Point>>("test laplacian_smoothing_for_unmarked_cells without move", answer_laplacian_smoothing_for_unmarked_cells<1>(), test_laplacian_smoothing_for_unmarked_cells<1>, [](std::vector<Point> a, std::vector<Point> b)->bool{return is_approximately_equal_Points(a, b);});
  test.test_approximate_value<std::vector<Point>>("test laplacian_smoothing_for_unmarked_cells with move", answer_laplacian_smoothing_for_unmarked_cells<2>(), test_laplacian_smoothing_for_unmarked_cells<2>, [](std::vector<Point> a, std::vector<Point> b)->bool{return is_approximately_equal_Points(a, b);});
  test.test("test laplacian_smoothing_for_unmarked_cells with Jacobians", test2_laplacian_smoothing_for_unmarked_cells<1>);
  test.test("test move_points_onto_mesh slanted plane", test_move_points_onto_mesh<1>);
  test.test("test move_points_onto_mesh 2 plane", test_move_points_onto_mesh<2>);
  test.test("test volumes_around_node", test_volumes_around_node<1>);
  test.test("test __get_signal", test_get_signal<1>);
  test.test("test __get_seven_non_manifold_templates", test_get_seven_non_manifold_templates<1>);
  test.test("test __get_non_manifold_template_list", test_get_non_manifold_template_list<1>);
  test.test("test __get_solution_to_non_manifold_templates_list", test_get_solution_to_non_manifold_templates_list<1>);
  test.test("test resolve_non_manifold_case with no problem", test_resolve_non_manifold_case<1>);
  test.test("test resolve_non_manifold_case with problem 1", test_resolve_non_manifold_case<2>);
  test.test("test resolve_non_manifold_case with problem 2", test_resolve_non_manifold_case<3>);
  // test.test("test resolve_non_manifold_case which should raise an error", test_resolve_non_manifold_case<4>);
  test.test("test detect_intersection_with_volume_fraction", test_detect_intersection_with_volume_fraction<1>);
  test.test("test move_points_onto_mesh_with_volume_fraction slanted plane", test_move_points_onto_mesh_with_volume_fraction<1>);
  test.test("test get_orthogonal_vectors", test_get_orthogonal_vectors<1>);
  test.test("test get_neighbors_list_for_smoothing slanted plane",test_get_neighbors_list_for_smoothing<1>);
  test.test("test get_neighbors_list_for_smoothing 2 planes and 4 cubes",test_get_neighbors_list_for_smoothing<2>);
  test.test("test surface_smoothing a plane", test_surface_smoothing<1>);
  test.test("test surface_smoothing a curved plane", test_surface_smoothing<2>);

  test.print_summary();
}