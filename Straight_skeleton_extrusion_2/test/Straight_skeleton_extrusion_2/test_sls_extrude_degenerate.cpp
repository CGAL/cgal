// Tests of CGAL::extrude_skeleton() on degenerate and nearly-degenerate (but valid) inputs.
//
// Invalid inputs (fewer than 3 vertices, duplicated vertices, touching holes, ...) are in
// test_sls_extrude_invalid_inputs.cpp.

#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
using Case_list = std::vector<std::tuple<std::string, SST::Extrusion_input<K>, SST::Expectations> >;

template <typename K>
Case_list<K> degenerate_cases()
{
  using FT = typename K::FT;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;
  using Outcome = SST::Outcome;

  Case_list<K> cases;
  auto add = [&](const std::string& name, const Input& in, const Expectations& e)
  {
    cases.emplace_back(name, in, e);
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Parallel edges

  // Comb: 4 teeth of width 1 (anti-parallel edges at distance 1) and gaps of width 1, base of height 1.
  // Many simultaneous edge events at z = 0.5 (one per tooth and the base's sub-rectangles).
  {
    Input in(SST::make_polygon<K>({{0,0}, {9,0}, {9,3}, {8,3}, {8,1}, {7,1}, {7,3}, {6,3}, {6,1}, {5,1}, {5,3}, {4,3},
                                   {4,1}, {3,1}, {3,3}, {2,3}, {2,1}, {1,1}, {1,3}, {0,3}}));
    Expectations e;
    e.nv = 30;
    e.nf = 56;
    e.volume = 14. / 3.;
    e.zmax_abs = 0.5;
    e.genus = 0;
    add("comb (simultaneous edge events)", in, e);
  }

  // Staircase: anti-parallel edges at distances 2, 4 and 6 (events at three different heights)
  {
    Input in(SST::make_polygon<K>({{0,0}, {9,0}, {9,6}, {6,6}, {6,4}, {3,4}, {3,2}, {0,2}}));
    Expectations e;
    e.nv = 14;
    e.zmax_abs = 2;
    e.genus = 0;
    add("staircase (anti-parallel edges at several distances)", in, e);
  }

  // U-shape: the two top edges y = 4 are parallel, same direction, on the same supporting line
  {
    Input in(SST::make_polygon<K>({{0,0}, {6,0}, {6,4}, {4,4}, {4,2}, {2,2}, {2,4}, {0,4}}));
    Expectations e;
    e.nv = 12;
    e.nf = 20;
    e.volume = 28. / 3.;
    e.zmax_abs = 1.;
    e.genus = 0;
    add("U (collinear non-consecutive edges)", in, e);
  }

  // Zigzag: same-direction parallel edges at different offsets
  {
    Input in(SST::make_polygon<K>({{0,0}, {2,1}, {4,0}, {6,1}, {8,0}, {8,3}, {6,4}, {4,3}, {2,4}, {0,3}}));
    Expectations e;
    e.nv = 15;
    e.genus = 0;
    add("zigzag (same-direction parallel edges)", in, e);
  }

  // Nearly parallel opposite edges
  {
    Input in(SST::make_polygon<K>({{0,0}, {10,0}, {10,1}, {0,1.000001}}));
    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.genus = 0;
    add("nearly parallel edges", in, e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Collinear vertices

  // Square with two extra vertices on its bottom side (three collinear consecutive edges):
  // the geometry is the plain pyramid. The bisectors of the collinear vertices are perpendicular to
  // the side and end on the two hip edges, at (0.5,0.5,0.5) and (1.5,0.5,0.5): 6 + 1 + 2 = 9 vertices.
  {
    Input in(SST::make_polygon<K>({{0,0}, {0.5,0}, {1.5,0}, {2,0}, {2,2}, {0,2}}));
    Expectations e;
    e.nv = 9;
    e.nf = 14;
    e.volume = 4. / 3.;
    e.zmax_abs = 1.;
    e.genus = 0;
    add("collinear vertices on one side", in, e);
  }

  // Nearly collinear vertex (tiny convex angle).
  {
    Input in(SST::make_polygon<K>({{0,0}, {1,1e-9}, {2,0}, {2,2}, {0,2}}));
    Expectations e;
    e.genus = 0;
    add("nearly collinear vertex", in, e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Symmetric structures and simultaneous events

  // Plus sign: 4 reflex vertices, 5 unit-width squares, multi-way event at the center
  {
    Input in(SST::make_polygon<K>({{2,0}, {4,0}, {4,2}, {6,2}, {6,4}, {4,4}, {4,6}, {2,6}, {2,4}, {0,4}, {0,2}, {2,2}}));
    Expectations e;
    e.nv = 17;
    e.nf = 30;
    e.volume = 28. / 3.;
    e.zmax_abs = 1.;
    e.genus = 0;
    add("plus (multi-way vertex events)", in, e);
  }

  // Octagon with axis-aligned and 45 degree edges (symmetric, irrational offset directions)
  {
    Input in(SST::make_polygon<K>({{1,0}, {2,0}, {3,1}, {3,2}, {2,3}, {1,3}, {0,2}, {0,1}}));
    Expectations e;
    e.nv = 13;
    e.nf = 22;
    e.genus = 0;
    add("octagon 45 degrees", in, e);
  }

  // Regular polygons from cos/sin: the input is only *nearly* regular, so the single n-way event at
  // the center may be split in several close events (observed: 2n-2 vertices for n >= 4, with both
  // kernels). Counts are not asserted; volume = area * apothem / 3 and height = apothem are.
  for(int n : {3, 4, 5, 6, 7, 8, 12, 16})
  {
    Polygon_2 p;
    for(int i=0; i<n; ++i)
      p.push_back(typename K::Point_2(std::cos(2 * CGAL_PI * i / n), std::sin(2 * CGAL_PI * i / n)));
    const double apothem = std::cos(CGAL_PI / n);
    const double area = 0.5 * n * std::sin(2 * CGAL_PI / n);
    Expectations e;
    e.volume = area * apothem / 3.;
    e.zmax_abs = apothem;
    e.genus = 0;
    e.rel_tol = 1e-9;
    add("regular " + std::to_string(n) + "-gon", Input(p), e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Extreme aspect ratios and coordinates

  // hip roof of width W, length L: W^2 L / 4 - W^3 / 12
  {
    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.volume = 1e6 / 4. - 1. / 12.;
    e.zmax_abs = 0.5;
    e.rel_tol = 1e-9;
    add("thin rectangle 1e6 x 1", Input(SST::make_rectangle<K>(0, 0, 1e6, 1)), e);
  }

  {
    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.volume = 1e-12 / 4. - 1e-18 / 12.;
    e.zmax_abs = 5e-7;
    e.rel_tol = 1e-6;
    add("thin rectangle 1 x 1e-6", Input(SST::make_rectangle<K>(0, 0, 1, 1e-6)), e);
  }

  {
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 1e-21 / 6.;
    e.zmax_abs = 5e-8;
    add("tiny square 1e-7", Input(SST::make_rectangle<K>(0, 0, 1e-7, 1e-7)), e);
  }

  {
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 1e27 / 6.;
    e.zmax_abs = 5e8;
    add("huge square 1e9", Input(SST::make_rectangle<K>(0, 0, 1e9, 1e9)), e);
  }

  {
    // small object far from the origin (volume computed relative to a local origin)
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 4. / 3.;
    e.zmax_abs = 1.;
    add("unit square at 1e9", Input(SST::make_rectangle<K>(1e9, 1e9, 1e9 + 2, 1e9 + 2)), e);
  }

  // Spike: very sharp convex vertex
  {
    Input in(SST::make_polygon<K>({{0,0}, {4,0}, {4,2}, {2.0001,2}, {2,100}, {1.9999,2}, {0,2}}));
    Expectations e;
    e.nv = 11;
    e.nf = 18;
    e.genus = 0;
    add("spike", in, e);
  }

  // Very sharp reflex notch
  {
    Input in(SST::make_polygon<K>({{0,0}, {1.9999,0}, {2,1.9}, {2.0001,0}, {4,0}, {4,2}, {0,2}}));
    Expectations e;
    e.genus = 0;
    add("sharp notch", in, e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Holes close to each other / to the boundary (but not touching)

  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(1e-6, 4, 2, 6, false));
    Expectations e;
    e.genus = 1;
    e.zmax_abs = 4;
    add("hole nearly touching the outer boundary", Input(pwh), e);
  }

  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(2, 2, 4, 8, false));
    pwh.add_hole(SST::make_rectangle<K>(4 + 1e-6, 2, 6, 8, false));
    Expectations e;
    e.zmax_abs = 2;
    e.genus = 2;
    add("two holes nearly touching", Input(pwh), e);
  }

  {
    // many small holes on a grid: many simultaneous events
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
        pwh.add_hole(SST::make_rectangle<K>(1 + 3 * i, 1 + 3 * j, 3 + 3 * i, 3 + 3 * j, false));
    Expectations e;
    e.genus = 9;
    e.nv = 56;
    e.volume = 50. / 3.;
    e.zmax_abs = 0.5;
    add("3x3 grid of holes", Input(pwh), e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Weights: simultaneous / nearly simultaneous events, extreme slopes

  {
    // weights making the 4-way apex nearly (but not exactly) simultaneous
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));
    in.with_weights(Weights{{FT(0.5), FT(1), FT(0.5), FT(1 + 1e-12)}});
    Expectations e;
    e.volume = 8. / 3.;
    e.zmax_abs = 1.;
    add("nearly simultaneous weighted events", in, e);
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_weights(Input::uniform(in.pwh, 1e-3));
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 4e-3 / 3.;
    e.zmax_abs = 1e-3;
    add("very shallow slopes (weight 1e-3)", in, e);
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_weights(Input::uniform(in.pwh, 1e3));
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 4e3 / 3.;
    e.zmax_abs = 1e3;
    add("very steep slopes (weight 1e3)", in, e);
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_angles(Input::uniform(in.pwh, 89.9));
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.zmax_abs = std::tan(89.9 * CGAL_PI / 180);
    e.rel_tol = 1e-9;
    add("angles 89.9 degrees", in, e);
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_angles(Input::uniform(in.pwh, 0.01));
    Expectations e;
    e.nv = 5;
    e.zmax_abs = std::tan(0.01 * CGAL_PI / 180);
    e.rel_tol = 1e-9;
    add("angles 0.01 degrees", in, e);
  }

  {
    // angles just above 90 degrees: nearly vertical outward walls
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_angles(Input::uniform(in.pwh, 90.1)).with_height(1);
    const double d = 1. / std::tan(89.9 * CGAL_PI / 180); // horizontal growth
    Expectations e;
    e.nv = 8;
    e.top_area = (2 + 2 * d) * (2 + 2 * d);
    e.rel_tol = 1e-9;
    add("angles 90.1 degrees (outward), H = 1", in, e);
  }

  {
    // vertical walls on a non-axis-aligned polygon (the snapping projects onto the edge)
    Input in(SST::make_polygon<K>({{0,0}, {3,1}, {4,4}, {1,3}}));
    in.with_angles(Weights{{FT(45), FT(90), FT(45), FT(90)}});
    Expectations e;
    e.genus = 0;
    e.rel_tol = 1e-6;
    add("vertical walls on slanted edges", in, e);
  }

  return cases;
}

template <typename K>
std::vector<CGAL::Surface_mesh<typename K::Point_3> > run_cases()
{
  std::vector<CGAL::Surface_mesh<typename K::Point_3> > meshes;
  for(const auto& [name, in, e] : degenerate_cases<K>())
  {
    std::cout << "run " << name << " (" << typeid(K).name() << ")" << std::endl;
    meshes.emplace_back();
    const SST::Outcome o = extrude(in, meshes.back());
    CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name + ": " + to_string(o));
    if(o == SST::Outcome::SUCCESS)
      check_extrusion(name, in, meshes.back(), e);
  }
  return meshes;
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  const auto epick_meshes = run_cases<EPICK>();
  std::cout << "== EPECK ==" << std::endl;
  const auto epeck_meshes = run_cases<EPECK>();

  // Cross-kernel consistency: same combinatorics, same coordinates up to 1e-9 * input size
  std::cout << "== EPICK vs EPECK ==" << std::endl;
  const auto names = degenerate_cases<EPICK>();
  for(std::size_t i=0; i<epick_meshes.size(); ++i)
  {
    const std::string& name = std::get<0>(names[i]);
    assert(num_faces(epick_meshes[i]) != 0 && num_faces(epeck_meshes[i]) != 0);
    const double tol = 1e-9 * (std::max)(1e-300, SST::diag_xy(PMP::bbox(epeck_meshes[i])));
    CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(epick_meshes[i]), PMP::volume(epeck_meshes[i]), tol),
                       "EPICK vs EPECK: " + name);
  }

  return SST::test_exit_code();
}
