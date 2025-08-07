#include "hexmeshing_sequential.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/draw_polyhedron.h>
#include <iostream>
#include <cassert>

using namespace CGAL::HexRefinement::TwoRefinement;

// A modifier creating a icosahedron with the incremental builder.
template <class HDS>
class Build_icosahedron : public CGAL::Modifier_base<HDS> {
public:
  Build_icosahedron() {}
  void operator()( HDS& hds) {
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    B.begin_surface( 12, 20, 0);

    const double phi2 = 1. + std::sqrt(5.);
    std::vector<Point> vertices = {
      Point(0, 2, phi2),
      Point(0, 2, -phi2),
      Point(0, -2, phi2),
      Point(0, -2, -phi2),
      Point(phi2, 0, 2),
      Point(-phi2, 0, 2),
      Point(phi2, 0, -2),
      Point(-phi2, 0, -2),
      Point(2, phi2, 0),
      Point(2, -phi2, 0),
      Point(-2, phi2, 0),
      Point(-2, -phi2, 0)
    };
    for (const auto& v : vertices) {
      B.add_vertex(v + Vector(5, 5, 5));
    }
    
    std::vector<std::vector<int>> faces = {
      {0, 4, 8}, {1, 7, 10}, {2, 5, 11}, {3, 6, 9},
      {1, 8, 6}, {0, 10, 5}, {2, 9, 4}, {3, 11, 7},
      {0, 2, 4}, {2, 0, 5}, {1, 3, 7}, {3, 1, 6},
      {4, 6, 8}, {6, 4, 9}, {5, 7, 11}, {7, 5, 10},
      {8, 10, 0}, {10, 8, 1}, {9, 11, 3}, {11, 9, 2}
    };
    
    for (const auto& f : faces) {
      B.begin_facet();
      for (int v : f) {
        B.add_vertex_to_facet(v);
      }
      B.end_facet();
    }
    B.end_surface();
  }
};


// A modifier creating a icosahedron with the incremental builder.
template <class HDS>
class Build_cube : public CGAL::Modifier_base<HDS> {
public:
  Build_cube() {}
  void operator()( HDS& hds) {
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    B.begin_surface(8, 12, 0);

    const double l = 2.8;
    std::vector<Point> vertices = {
      Point(l, l, l),
      Point(l, l, -l),
      Point(l, -l, l),
      Point(l, -l, -l),
      Point(-l, l, l),
      Point(-l, l, -l),
      Point(-l, -l, l),
      Point(-l, -l, -l)
    };
    for (const auto& v : vertices) {
      B.add_vertex(v + Vector(5, 5, 5));
    }
    
    std::vector<std::vector<int>> faces = {
      {0, 1, 3}, {0, 3, 2}, {4, 6, 7}, {4, 7, 5}, {0, 2, 6}, {0, 6, 4}, {1, 5, 7}, {1, 7, 3}, {0, 4, 5}, {0, 5, 1}, {7, 6, 2}, {7, 2, 3}
    };
    
    for (const auto& f : faces) {
      B.begin_facet();
      for (int v : f) {
        B.add_vertex_to_facet(v);
      }
      B.end_facet();
    }
    B.end_surface();
  }
};

using HalfedgeDS = Polyhedron::HalfedgeDS;

template<int FACE=4>
Polyhedron create_polyhedron() {
  Polyhedron poly;
  Build_icosahedron<HalfedgeDS> icosahedron;
  Build_cube<HalfedgeDS> cube;
  const double c = 5.6;
  const double d = 2.4;
  switch(FACE) {
    case 20:
      poly.delegate(icosahedron);
      break;
    case 6:
      poly.delegate(cube);
      break;
    case 4:
      poly.make_tetrahedron(
        Point(c+d, c+d, c+d),
        Point(c-d, c+d, c-d),
        Point(c+d, c-d, c-d),
        Point(c-d, c-d, c+d)
      );
      break;
    default:
      std::cout << "NO POLYHEDRON GENERATED" << std::endl;
  }
  return poly;
}

Tree get_surface_aabb(Polyhedron& poly) {
  // Triangulate before AABB
  Tree tree;
  // Compute AABB tree
  tree.insert(faces(poly).first, faces(poly).second, poly);
  tree.accelerate_distance_queries();

  return tree;
}

// テスト用のrefinementされたメッシュ作成
LCC create_refined_test_mesh(int level = 1) {
  std::cout << "Creating refined test mesh..." << std::endl;
  
  Grid grid = Grid::make_cube(Point(0, 0, 0), 1.0, 10);

  Polyhedron poly = create_polyhedron();
  Tree aabb = get_surface_aabb(poly);
  
  return CGAL::HexRefinement::two_refinement(grid,
                        is_volume_intersecting_poly(aabb),
                        is_inner_point(aabb),
                        level,
                        false);
}

template<int FACE=6>
LCC create_refined_test_mesh_with_volume_fraction() {
  std::cout << "Creating refined test mesh..." << std::endl;
  
  Grid grid = Grid::make_cube(Point(0, 0, 0), 1.0, 13);

  LCC lcc;
  generate_grid(lcc, grid);
  
  const double l0 = 0.7;

  set_centroids(lcc);
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto it = volumes.begin(); it != volumes.end(); it++) {
    auto &vol_attr = lcc.attribute<3>(it)->info();
    Point c = vol_attr.centroid;
    double f = 1.;
    std::array<double, 3> xyz = {c.x(), c.y(), c.z()};
    if constexpr(FACE == 6) {
      for(auto x:xyz) {
        if(3 < x and x < 10) f *= 1.;
        else if(x < 2 or 11 < x) f *= 0.;
        else f *= l0;
      }
    }
    else if constexpr(FACE == 1) {
      const double maku = 5.5;
      double s = 0;
      for(auto x:xyz) {
        s += int(x);
      }
      if(maku < s) f = 0.;
      else if(maku < s+1) {
        double t = maku-s;
        f = t*t*t/6.;
      }
      else if(maku < s+2) {
        double t = maku-s, t1 = t-1.;
        f = t*t*t/6. - t1*t1*t1/2.;
      }
      else if(maku < s+3) {
        double t = s+3-maku;
        f = 1. - t*t*t/6.;
      }
    }
    vol_attr.fraction = f;
  }
  return lcc;
}

// メッシュへの点移動のテスト
void test_move_points_onto_mesh() {
  std::cout << "=== Testing Move Points Onto Mesh ===" << std::endl;
  
  LCC lcc = create_refined_test_mesh(0);
  
  Polyhedron poly = create_polyhedron();
  Tree aabb = get_surface_aabb(poly);
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, true);
  
  size_type move_mark = lcc.get_new_mark();
  
  // 移動前の座標を記録
  std::vector<Point> original_points;
  auto vertices = lcc.one_dart_per_cell<0>();
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    original_points.push_back(lcc.point(it));
  }
  
  auto func = detect_intersection(aabb, poly);
  // 正二十面体との交差検出関数
  // set_dual_edges(lcc);
  // auto voxes = lcc.one_dart_per_cell<3>();
  // int cnt = 0;
  // for(auto it = voxes.begin(); it != voxes.end(); it++, cnt++) {
  //   auto &face_attr = lcc.attribute<2>(it)->info();
  //   if(func(lcc, it)) {
  //     std::cout << "dual_edge: " << face_attr.dual_edge << std::endl;
  //     std::cout << "face_attr.intersection: " << face_attr.intersection << std::endl;
  //     std::cout << "face_attr.normal: " << face_attr.normal << std::endl;
  //   }
  // }
  
  // 点移動を実行
  move_points_onto_mesh(lcc, move_mark, func);
  
  // 移動後の座標を確認
  int moved_count = 0;
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    if (lcc.is_marked(it, move_mark)) {
      auto id = get_or_create_attr<0>(lcc, it)->id;
      Point moved_point = lcc.point(it);
      Point original_point = original_points[id];
      
      std::cout << "Vertex " << id << " moved from (" 
            << original_point.x() << ", " << original_point.y() << ", " << original_point.z() 
            << ") to (" << moved_point.x() << ", " << moved_point.y() << ", " << moved_point.z() << ")" << std::endl;
      
      // 座標が有効な範囲内にあることを確認
      assert(moved_point.x() >= 0 && moved_point.x() <= 10);
      assert(moved_point.y() >= 0 && moved_point.y() <= 10);
      assert(moved_point.z() >= 0 && moved_point.z() <= 10);
      
      moved_count++;
    }
  }
  
  std::cout << "Move points onto mesh test passed for " << moved_count << " vertices!" << std::endl;
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  auto trimming_func = is_inner_centroid(aabb);
  // auto trimming_func = is_volume_intersecting_poly(aabb);
  trim_excedent_volumes(lcc, trimming_func);
  CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  
  // surface_smoothing(lcc, move_mark);
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);

  lcc.free_mark(move_mark);
}

// 最初から volume fraction で与えられるパターンに対応したい
void test_move_points_onto_mesh_with_volume_fraction() {
  std::cout << "=== Testing Move Points Onto Mesh With Volume Fraction ===" << std::endl;

  const int refinement_level = 1;

  Polyhedron poly = create_polyhedron();
  Tree aabb = get_surface_aabb(poly);
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, true);
  
  auto cellIdentifier = is_volume_intersecting_poly(aabb);
  auto decideFunc = is_inner_point(aabb);

  // LCC lcc = create_refined_test_mesh(refinement_level);
  // set_fraction(lcc, 1./(1<<refinement_level), cellIdentifier, decideFunc);
  LCC lcc = create_refined_test_mesh_with_volume_fraction();

  size_type half_mark = lcc.get_new_mark();
  set_centroids(lcc);
  auto volumes = lcc.one_dart_per_cell<3>();
  for(auto volume = volumes.begin(); volume != volumes.end(); volume++) {
    // if(0.1 < lcc.attribute<3>(volume)->info().fraction and lcc.attribute<3>(volume)->info().fraction < 0.9) {
    //   std::cout << lcc.attribute<3>(volume)->info().fraction << ' ' << lcc.attribute<3>(volume)->info().centroid << std::endl;
    // }
    if(lcc.attribute<3>(volume)->info().centroid.z() < 5.0)
      lcc.mark_cell<3>(volume, half_mark);
  }

  // 移動前の座標を記録
  std::vector<Point> original_points;
  auto vertices = lcc.one_dart_per_cell<0>();
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    original_points.push_back(lcc.point(it));
  }
  
  size_type move_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();
  move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);
  
  // 移動後の座標を確認
  int moved_count = 0;
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    if (lcc.is_marked(it, move_mark)) {
      auto id = get_or_create_attr<0>(lcc, it)->id;
      Point moved_point = lcc.point(it);
      Point original_point = original_points[id];
      
      std::cout << "Vertex " << id << " moved from (" 
            << original_point.x() << ", " << original_point.y() << ", " << original_point.z() 
            << ") to (" << moved_point.x() << ", " << moved_point.y() << ", " << moved_point.z() << ")" << std::endl;
      
      // 座標が有効な範囲内にあることを確認
      // assert(moved_point.x() >= 0 && moved_point.x() <= 10);
      // assert(moved_point.y() >= 0 && moved_point.y() <= 10);
      // assert(moved_point.z() >= 0 && moved_point.z() <= 10);
      
      moved_count++;
    }
  }
  
  std::cout << "Move points onto mesh test passed for " << moved_count << " vertices!" << std::endl;
  // CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  auto trimming_func = is_marked_volume(inner_mark);
  // auto trimming_func = is_volume_intersecting_poly(aabb);
  trim_excedent_volumes(lcc, trimming_func);
  std::cout << "Rendering the result of move_points_onto_mesh_with_volume_fraction" << std::endl;
  CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);
  
  surface_smoothing(lcc, move_mark, inner_mark);
  std::cout << "Rendering the result of surface_smoothing" << std::endl;
  CGAL::HexRefinement::render_two_refinement_result(lcc, aabb, false);

  CGAL::HexRefinement::render_two_refinement_result_with_mark(lcc, half_mark);
  volume_smoothing(lcc, move_mark);
  std::cout << "Rendering the result of volume_smoothing" << std::endl;
  CGAL::HexRefinement::render_two_refinement_result_with_mark(lcc, half_mark);

  lcc.free_mark(move_mark);
  lcc.free_mark(inner_mark);
  lcc.free_mark(half_mark);
}

int main() {
  // test_move_points_onto_mesh();
  test_move_points_onto_mesh_with_volume_fraction();
}