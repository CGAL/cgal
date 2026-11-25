// Tests of CGAL::extrude_skeleton() on degenerate and nearly-degenerate (but valid) inputs.
//
// Invalid inputs (fewer than 3 vertices, duplicated vertices, touching holes, ...) are in
// test_sls_extrude_invalid_inputs.cpp.

#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
using Case_list = std::vector<std::tuple<std::string, SST::Extrusion_input<K>, SST::Expectations> >;

// Weighted polygon with holes generated as in Straight_skeleton_2/test/Straight_skeleton_2/test_sls_weighted_polygons_with_holes.cpp:
// a long rectangle with fast vertical sides, and 2 * `hole_n` random holes with slow edges. The skeleton faces
// of the vertical sides surround the faces of the holes, so the builder adds artificial bisectors.
std::pair<CGAL::Polygon_with_holes_2<EPICK>, std::vector<std::vector<double> > >
generate_holed_rectangle(const int hole_n, const int hole_nv, const unsigned int seed)
{
  using Point_2 = EPICK::Point_2;
  using Polygon_2 = CGAL::Polygon_2<EPICK>;

  CGAL::Random rnd(seed);

  std::vector<Point_2> ob = { Point_2(-hole_n-1, -0.5), Point_2(hole_n+1, -0.5), Point_2(hole_n+1, 0.5), Point_2(-hole_n-1, 0.5) };
  CGAL::Polygon_with_holes_2<EPICK> pwh(Polygon_2(ob.begin(), ob.end()));

  std::vector<std::vector<double> > weights(1);
  weights[0].push_back(rnd.get_double(0.05, 0.5));
  weights[0].push_back(rnd.get_double(1, 10));
  if(rnd.get_int(0, 2) == 1)
    weights[0].push_back(rnd.get_double(0.05, 0.5));
  else
    weights[0].push_back(weights[0][0]);
  weights[0].push_back(rnd.get_double(1, 10));

  for(int hole_i=0; hole_i<hole_n; ++hole_i)
  {
    for(const int side : {1, -1})
    {
      Polygon_2 poly;
      CGAL::random_polygon_2(hole_nv, std::back_inserter(poly), CGAL::Random_points_in_square_2<Point_2>(0.25, rnd));
      std::vector<Point_2> hole_pts = poly.container();
      const double y_nudge = rnd.get_double(-0.25, 0.25);
      for(Point_2& pt : hole_pts)
        pt = Point_2(pt.x() + side * (hole_i + 1), pt.y() + y_nudge);
      // holes must be clockwise; same as Polygon_2::reverse_orientation(), which keeps the first vertex
      if(CGAL::polygon_area_2(hole_pts.begin(), hole_pts.end(), EPICK()) > 0)
        std::reverse(std::next(hole_pts.begin()), hole_pts.end());
      pwh.add_hole(Polygon_2(hole_pts.begin(), hole_pts.end()));
      weights.emplace_back(hole_pts.size(), 10.);
    }
  }

  return { pwh, weights };
}

// Maximal nesting of artificial bisectors in a skeleton face: 1 if the artificial bisectors link holes to the border
// of the face, 2 if a hole is linked to the border through another hole, etc.
// Artificial bisectors are the only halfedges with the same face on both sides, and the boundary cycle of a face
// goes along each artificial bisector twice, first towards the hole and then back.
template <typename Ss>
int artificial_bisector_nesting(const Ss& ss)
{
  int max_depth = 0;
  for(auto f = ss.faces_begin(); f != ss.faces_end(); ++f)
  {
    std::vector<int> stack; // the vertex where we left each enclosing loop
    int last = -1;
    auto h = f->halfedge();
    do
    {
      const int vid = int(h->vertex()->id());
      if(h->face() == h->opposite()->face())
      {
        if(!stack.empty() && stack.back() == vid)
          stack.pop_back();
        else
          stack.push_back(last);
        max_depth = (std::max)(max_depth, int(stack.size()));
      }
      last = vid;
      h = h->next();
    }
    while(h != f->halfedge());
  }
  return max_depth;
}

template <typename K>
Case_list<K> degenerate_cases()
{
  using FT = typename K::FT;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;

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

  // The faces of the left and right edges meet along a horizontal ridge (its vertices have the same time),
  // and the height is the time of the ridge, as computed with doubles: the top face is a sliver of (almost)
  // zero width around the ridge, which can be inverted by rounding.
  {
    const std::vector<std::pair<std::vector<double>, double> > ws_and_hs = {
      { { 0.48985014453064651, 5.7434052755124867, 0.27788849917706104, 1.7344466173090041 }, 0.35460953437980813 },
      { { 0.11170436046086252, 2.1415360951796174, 0.20877599539235236, 9.0026620849967003 }, 0.14553896124331675 } };
    for(const auto& [ws, H] : ws_and_hs)
    {
      Input in(SST::make_rectangle<K>(-1, -0.5, 1, 0.5));
      in.with_weights(Weights{{FT(ws[0]), FT(ws[1]), FT(ws[2]), FT(ws[3])}}).with_height(H);
      Expectations e;
      e.genus = 0;
      e.top_area = 0; // the sign of the (tiny) area of the top face is not reliable with EPICK
      e.rel_tol = 1e-6;
      add("height at the time of a horizontal ridge, H = " + std::to_string(H), in, e);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Non-simply connected skeleton faces

  // The fast right edge sweeps around the hole: its face surrounds the faces of the hole. The builder links them
  // with an artificial bisector, from (5.56, 1.06, t = 0.44) to (0.91, 1.06, t = 0.91), which the (single) boundary
  // of the face goes along in both directions. H = 0.3 is below the artificial bisector, H = 0.45 cuts it.
  for(double H : {0., 0.3, 0.45, -0.6})
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 4));
    pwh.add_hole(SST::make_rectangle<K>(6, 1.5, 7, 2.5, false /*cw*/));
    Input in(pwh);
    in.with_weights(Weights{{FT(1), FT(1), FT(0.1), FT(1)}, {FT(1), FT(1), FT(1), FT(1)}});
    if(H != 0)
      in.with_height(H);
    Expectations e;
    e.genus = 1;
    e.rel_tol = 1e-6;
    add("annular skeleton face, H = " + std::to_string(H), in, e);
  }

  // An artificial bisector that ends at an actual skeleton vertex, or very close to one. The fast right edge sweeps
  // around the hole (near the right edge), and the artificial bisector is shot from the leftmost vertex of the faces
  // of the hole along y = 2, until it hits the border of the face of the right edge:
  //  - w = 0.3: at the vertex (3.33, 2, t = 2) where the faces of the right, top, and bottom edges meet;
  //  - w = 0.25: at the vertex (2, 2, t = 2) where these faces and the face of the left edge meet;
  //  - w = 0.25 - 1e-12: on the bisector of the left and right edges, at about 6e-12 from (2, 2, t = 2).
  // H = 2 is the time of these vertices.
  for(const double w : {0.3, 0.25, 0.25 - 1e-12})
  {
    // Reference skeleton, to check that the configuration is as described
    {
      CGAL::Polygon_with_holes_2<EPICK> i_pwh(SST::make_rectangle<EPICK>(0, 0, 10, 4));
      i_pwh.add_hole(SST::make_polygon<EPICK>({{8.5, 2}, {9, 2.5}, {9.5, 2}, {9, 1.5}}));
      std::vector<std::vector<double> > i_weights = {{1, 1, w, 1}, {1, 1, 1, 1}};
      const auto ss = CGAL::create_interior_weighted_straight_skeleton_2(i_pwh, i_weights, EPICK());
      CGAL_SLS_CHECK_MSG(ss != nullptr, "w = " << w << ": no skeleton");
      if(!ss)
        continue;

      int artificial_bisectors = 0, to_actual_vertex = 0;
      for(auto h = ss->halfedges_begin(); h != ss->halfedges_end(); ++h)
      {
        if(h->face() != h->opposite()->face() || h->opposite() < h)
          continue;
        ++artificial_bisectors;
        // the end that is not the extreme vertex of the hole is the farthest from the contour edge
        auto v = (h->vertex()->time() > h->opposite()->vertex()->time()) ? h->vertex() : h->opposite()->vertex();
        if(v->trisegment()->e0().id() != v->trisegment()->e1().id()) // not the marker of artificial vertices
          ++to_actual_vertex;
      }
      CGAL_SLS_CHECK_MSG(artificial_bisectors == 1 && to_actual_vertex == ((w == 0.25 - 1e-12) ? 0 : 1),
                         "w = " << w << ": " << artificial_bisectors << " artificial bisector(s), "
                                << to_actual_vertex << " to an actual vertex");
    }

    for(const double H : {0., 0.5, 1., 2.})
    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 4));
      pwh.add_hole(SST::make_polygon<K>({{8.5, 2}, {9, 2.5}, {9.5, 2}, {9, 1.5}})); // clockwise
      Input in(pwh);
      in.with_weights(Weights{{FT(1), FT(1), FT(w), FT(1)}, {FT(1), FT(1), FT(1), FT(1)}});
      if(H != 0)
        in.with_height(H);
      Expectations e;
      e.genus = 1;
      e.rel_tol = 1e-6;
      add("artificial bisector to an actual vertex (w = " + std::to_string(w) + "), H = " + std::to_string(H), in, e);
    }
  }

  // Several artificial bisectors in the same face: chained through holes (hole -> hole -> ... -> border of the face),
  // and several of them incident to the same hole. The face is then cut into an outer border and nested holes.
  //  - seed 33: nesting depth 3, and holes incident to three artificial bisectors
  //  - seed 47: nesting depth 6
  for(const unsigned int seed : {33u, 47u})
  {
    const auto [i_pwh, i_weights] = generate_holed_rectangle(4, 4, seed);

    // Reference skeleton, to check that the configuration is as described, and to get the maximal time
    std::vector<std::vector<EPICK::FT> > i_weights_ft(i_weights.begin(), i_weights.end());
    const auto ss = CGAL::create_interior_weighted_straight_skeleton_2(i_pwh, i_weights_ft, EPICK());
    CGAL_SLS_CHECK_MSG(ss != nullptr, "seed " << seed << ": no skeleton");
    if(!ss)
      continue;
    const int depth = artificial_bisector_nesting(*ss);
    CGAL_SLS_CHECK_MSG(depth >= 3, "seed " << seed << ": nesting of artificial bisectors is only " << depth);

    double t_max = 0;
    for(auto v = ss->vertices_begin(); v != ss->vertices_end(); ++v)
      t_max = (std::max)(t_max, CGAL::to_double(v->time()));

    CGAL::Cartesian_converter<EPICK, K> to_k;
    Polygon_2 outer;
    for(const auto& p : i_pwh.outer_boundary())
      outer.push_back(to_k(p));
    Polygon_with_holes_2 pwh(outer);
    for(const auto& h : i_pwh.holes())
    {
      Polygon_2 hole;
      for(const auto& p : h)
        hole.push_back(to_k(p));
      pwh.add_hole(hole);
    }

    Weights weights;
    for(const auto& ws : i_weights)
      weights.emplace_back(ws.begin(), ws.end());

    // H = 0: no height. Heights are generic: they do not coincide with the time of a skeleton vertex
    for(const double f : {0., 0.25, 0.5, 0.75})
    {
      Input in(pwh);
      in.with_weights(weights);
      if(f != 0)
        in.with_height(f * t_max);
      Expectations e;
      e.genus = int(pwh.number_of_holes());
      e.rel_tol = 1e-6;
      add("nested artificial bisectors (seed " + std::to_string(seed) + "), " +
          ((f == 0) ? std::string("no height") : "H = " + std::to_string(f) + " * t_max"), in, e);
    }
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
