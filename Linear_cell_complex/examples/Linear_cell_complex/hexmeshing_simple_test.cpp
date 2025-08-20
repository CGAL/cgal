#include <CGAL/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <CGAL/Hexmeshing_mesh_data_for_hexmeshing.h>
#include <CGAL/Hexmeshing_render_results.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/draw_polyhedron.h>
#include <iostream>
#include <cassert>


// A modifier creating a icosahedron with the incremental builder.
template <class HDS>
class Build_icosahedron : public CGAL::Modifier_base<HDS> {
public:
  Build_icosahedron() {}
  void operator()( HDS& hds) {
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

    using Vertex = typename HDS::Vertex;
    using Point = typename Vertex::Point;

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
      B.add_vertex(v + CGAL::Hexmeshing::Vector(5, 5, 5));
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

    using Vertex = typename HDS::Vertex;
    using Point = typename Vertex::Point;

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
      B.add_vertex(v + CGAL::Hexmeshing::Vector(5, 5, 5));
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

using namespace CGAL::Hexmeshing;
using HalfedgeDS = Polyhedron::HalfedgeDS;

template <int FACE=4>
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

// テスト用のrefinementされたメッシュ作成
template <int FACE>
LCC create_refined_test_mesh(int level = 1) {
  std::cout << "Creating refined test mesh..." << std::endl;
  
  Grid grid = Grid::make_cube(Point(0, 0, 0), 1.0, 10);

  Polyhedron poly = create_polyhedron<FACE>();
  CGAL::MeshDataForHexmeshing mesh(poly, grid);

  CGAL::HexMeshingData hdata;

  hdata.two_refinement(
    mesh,
    level,
    false
  );
  
  return hdata.lcc;
}

template<int FACE=6>
LCC create_refined_test_mesh_with_volume_fraction() {
  std::cout << "Creating refined test mesh..." << std::endl;
  
  Grid grid = Grid::make_cube(Point(0, 0, 0), 1.0, 13);

  LCC lcc;
  Grid::generate_grid(lcc, grid);
  
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

template <int FACE, bool start_with_volume_fraction>
void render_meshes_at_each_phase() {
  CGAL::HexMeshingData hdata;
  LCC& lcc = hdata.lcc;

  if constexpr(start_with_volume_fraction) {
    lcc = create_refined_test_mesh_with_volume_fraction<FACE>();
  } else {
    Polyhedron poly = create_polyhedron<FACE>();
    CGAL::MeshDataForHexmeshing mesh(poly);
    auto cellIdentifier = is_volume_intersecting_poly(*mesh.get_tree_pointer());
    auto decideFunc = is_inner_point(*mesh.get_tree_pointer());
    const int refinement_level = 1;
    lcc = create_refined_test_mesh<FACE>(refinement_level);
    set_fraction(lcc, 1./(1<<refinement_level), cellIdentifier, decideFunc);
  }

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
  
  size_type move_mark = lcc.get_new_mark();
  size_type inner_mark = lcc.get_new_mark();
  move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);
  
  auto trimming_func = is_marked_volume(inner_mark);
  // auto trimming_func = is_volume_intersecting_poly(aabb);
  trim_excedent_volumes(lcc, trimming_func);
  std::cout << "Rendering the result of move_points_onto_mesh_with_volume_fraction" << std::endl;
  render_two_refinement_result(hdata);
  
  surface_smoothing(lcc, move_mark, inner_mark);
  std::cout << "Rendering the result of surface_smoothing" << std::endl;
  render_two_refinement_result(hdata);

  render_two_refinement_result_with_mark(hdata, half_mark);
  volume_smoothing(lcc, move_mark);
  std::cout << "Rendering the result of volume_smoothing" << std::endl;
  render_two_refinement_result_with_mark(hdata, half_mark);

  lcc.free_mark(move_mark);
  lcc.free_mark(inner_mark);
  lcc.free_mark(half_mark);
}

int main() {
  render_meshes_at_each_phase<4, false>();
  render_meshes_at_each_phase<6, false>();
  render_meshes_at_each_phase<20, false>();
  render_meshes_at_each_phase<1, true>();
  render_meshes_at_each_phase<6, true>();
}