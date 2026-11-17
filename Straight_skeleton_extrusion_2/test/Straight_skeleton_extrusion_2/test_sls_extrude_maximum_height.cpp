#include "sls_extrusion_test_helpers.h"

#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/extrude_skeleton.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <set>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
struct Max_height_tests
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;
  using Outcome = SST::Outcome;

  static Mesh run(const std::string& name,
                  const Input& in,
                  const Expectations& e = Expectations())
  {
    std::cout << "run " << name << std::endl;

    Mesh m;
    const Outcome o = extrude(in, m);
    CGAL_SLS_CHECK_MSG(o == Outcome::SUCCESS, name + ": " + to_string(o));
    if(o == Outcome::SUCCESS)
      check_extrusion(name, in, m, e);
    return m;
  }

  // Area of the (unweighted) interior offset of `pwh` at distance d, from the 2D offset API
  static double offset_area(const Polygon_with_holes_2& pwh, const FT d)
  {
    auto offsets = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(d, pwh, K(), K());
    double a = 0;
    for(const auto& p : offsets)
      a += SST::polygon_area(*p);
    return a;
  }

  // Square [0,2]^2, uniform weights. Skeleton height = 1.
  // For H < 1: frustum, volume 4H - 4H^2 + 4H^3/3, top cap (2-2H)^2, 8 vertices / 12 faces.
  // For H >= 1: identical to the unbounded pyramid.
  static void test_square_heights()
  {
    auto frustum_volume = [](const double H) {
      return (H >= 1) ? 4. / 3. : 4 * H - 4 * H * H + 4 * H * H * H / 3;
    };

    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);

    for(double H : {1e-12, 1e-9, 1e-6, 0.25, 0.5, 0.75, 1. - 1e-6, 1. - 1e-9})
    {
      Input in(sq);
      in.with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.genus = 0;
      e.volume = frustum_volume(H);
      e.zmax_abs = H;
      e.rel_tol = 1e-6;
      e.top_area = (2 - 2 * H) * (2 - 2 * H);
      run("square, H = " + std::to_string(H), in, e);
    }

    // equal to the skeleton height: no top cap, same as unbounded
    // (and above it, including the 'huge' values
    for(double H : {1., 1. + 1e-9, 2., 1e6, 1e30})
    {
      Input in(sq);
      in.with_height(H);
      Expectations e;
      e.nv = 5;
      e.nf = 6;
      e.genus = 0;
      e.volume = 4. / 3.;
      e.zmax_abs = 1.;
      run("square, H = " + std::to_string(H), in, e);
    }
  }

  // Rectangle [0,4]x[0,2]: ridge at z = 1.
  static void test_rectangle_ridge_coincidence()
  {
    const Polygon_2 r = SST::make_rectangle<K>(0, 0, 4, 2);
    {
      Input in(r);
      in.with_height(1);
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = 10. / 3.;
      e.top_area = 0.;
      run("rectangle, H = ridge height", in, e);
    }

    for(double eps : {1e-12, 1e-9, 1e-3})
    {
      const double H = 1 - eps;
      Input in(r);
      in.with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.zmax_abs = H;
      e.top_area = (4 - 2 * H) * (2 - 2 * H);
      e.volume = 8 * H - 6 * H * H + 4 * H * H * H / 3; // integral of (4-2t)(2-2t)
      e.rel_tol = 1e-6;
      run("rectangle, H = ridge - " + std::to_string(eps), in, e);
    }

    {
      Input in(r);
      in.with_height(1 + 1e-9);
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = 10. / 3.;
      run("rectangle, H = ridge + 1e-9", in, e);
    }
  }

  // Square ring [0,10]^2 \ [4,6]^2, consumed at t = 2. A(t) = 96 - 48 t.
  static void test_ring_heights()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));
    auto vol = [](double H) {
      H = (std::min)(H, 2.);
      return 96 * H - 24 * H * H;
    };

    for(double H : {0.5, 1., 1.5, 2. - 1e-9})
    {
      Input in(pwh);
      in.with_height(H);
      Expectations e;
      e.genus = 1;
      e.nv = 16;
      e.nf = 32;
      e.volume = vol(H);
      e.top_area = 96 - 48 * H;
      e.rel_tol = 1e-6;
      run("ring, H = " + std::to_string(H), in, e);
    }

    for(double H : {2., 3., 100.}) // consumed exactly at / before the cut
    {
      Input in(pwh);
      in.with_height(H);
      Expectations e;
      e.genus = 1;
      e.nv = 12;
      e.nf = 24;
      e.volume = 96.;
      run("ring, H = " + std::to_string(H), in, e);
    }
  }

  // Holes consumed at different heights (1, 1.25, 1.5): cut before, between, at and after.
  static void test_multiple_holes_heights()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(1, 1, 3, 3, false));
    pwh.add_hole(SST::make_rectangle<K>(6, 1, 9, 4, false));
    pwh.add_hole(SST::make_rectangle<K>(2, 6, 8, 8, false));

    for(double H : {0.25, 0.5, 0.75, 0.999, 1., 1.001, 1.25, 1.4, 1.5, 2.})
    {
      Input in(pwh);
      in.with_height(H);
      Expectations e;
      e.genus = 3;
      e.top_area = offset_area(pwh, FT(H));
      e.rel_tol = 1e-6;
      Mesh m = run("three holes, H = " + std::to_string(H), in, e);
      Mesh m2 = run("three holes, unbounded", Input(pwh));
      if(num_faces(m) > 0)
        CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(m), SST::clipped_volume<K>(m2, H), 1e-6),
                           "three holes, H = " << H);
    }
  }

  // Generic sweep: for each distinct height of an interior skeleton vertex (event height), cut at
  // that height and slightly below/above, and compare with the clipped unbounded extrusion.
  static void sweep_event_heights(const std::string& name,
                                  const Polygon_with_holes_2& pwh)
  {
    Mesh full = run(name + ", unbounded", Input(pwh));
    if(num_faces(full) == 0)
      return;

    std::set<FT> event_heights;
    for(auto v : vertices(full))
      if(CGAL::is_positive(full.point(v).z()))
        event_heights.insert(full.point(v).z());

    std::set<FT> heights; // sorted, so that the volume must be non-decreasing
    for(const FT& ze : event_heights)
    {
      const double zd = CGAL::to_double(ze);
      heights.insert(FT(zd * (1 - 1e-6)));
      heights.insert(ze);
      heights.insert(FT(zd * (1 + 1e-6)));
    }

    double prev_vol = 0;
    for(const FT& H : heights)
    {
      const std::string tname = name + ", H = " + std::to_string(CGAL::to_double(H)) +
                                  (event_heights.count(H) ? " (event)" : "");
      Input in(pwh);
      in.with_height(H);
      Expectations e;
      e.rel_tol = 1e-6;
      e.genus = int(pwh.number_of_holes());
      e.top_area = offset_area(pwh, H);
      Mesh m = run(tname, in, e);
      if(num_faces(m) == 0 || !CGAL::is_closed(m))
        continue;

      const FT v = PMP::volume(m);
      CGAL_SLS_CHECK_MSG(SST::approx_equal(v, SST::clipped_volume<K>(full, CGAL::to_double(H)), 1e-6, 1e-9), tname);
      CGAL_SLS_CHECK_MSG(v >= prev_vol * (1 - 1e-9), tname << ": volume not monotonic in H");
      prev_vol = CGAL::to_double(v);
    }
    CGAL_SLS_CHECK_MSG(SST::approx_equal(prev_vol, PMP::volume(full), 1e-5), name + ": last cut != full volume");
  }

  static void test_event_height_sweeps()
  {
    // two reflex vertices moving towards each other (vertex event / pinch of the offset curve).
    sweep_event_heights("pinch", Polygon_with_holes_2(SST::make_polygon<K>({{0,0}, {1,0}, {2,1}, {3,0}, {4,0}, {4,4},
                                                                            {3,4}, {2,3}, {1,4}, {0,4}})));

    // L with legs of different widths (events at z = 1 and z = 2)
    sweep_event_heights("uneven L", Polygon_with_holes_2(SST::make_polygon<K>({{0,0}, {6,0}, {6,2}, {4,2}, {4,6}, {0,6}})));

    // plus sign: 4 reflex vertices, multi-way event at the center
    sweep_event_heights("plus", Polygon_with_holes_2(SST::make_polygon<K>({{2,0}, {4,0}, {4,2}, {6,2}, {6,4}, {4,4},
                                                                           {4,6}, {2,6}, {2,4}, {0,4}, {0,2}, {2,2}})));

    // has a non-trivial ridge at z=1
    sweep_event_heights("ridge", Polygon_with_holes_2(SST::make_polygon<K>({{0,2}, {2,2}, {2,0}, {6,0}, {6,4}, {4,4},
                                                                            {4,6}, {0,6}})));

    // generic non-convex polygon (non axis-aligned edges).
    sweep_event_heights("star", Polygon_with_holes_2(SST::make_polygon<K>({{0,3}, {-1,1}, {-3,1}, {-1.5,-0.5}, {-2,-3},
                                                                           {0,-1.5}, {2,-3}, {1.5,-0.5}, {3,1}, {1,1}})));

    // polygon with an off-center hole (hole consumed on one side first)
    Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0,0}, {8,0}, {8,6}, {0,6}}));
    pwh.add_hole(SST::make_polygon<K>({{1,1}, {1,3}, {4,3}, {4,1}}));
    sweep_event_heights("off-center hole", pwh);
  }

  // Maximum height combined with weights, angles and vertical walls
  static void test_height_with_weights()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    {
      // one vertical wall (edge y = 0): region at t is [t,2-t] x [0,2-t]
      Input in(sq);
      in.with_weights(Weights{{FT(1), FT(0), FT(1), FT(1)}}).with_height(0.5);
      Expectations e;
      e.volume = 4. / 3.;
      e.top_area = 1. * 1.5;
      e.rel_tol = 1e-6;
      run("one vertical wall, H = 0.5", in, e);
    }

    {
      // all walls vertical (a prism)
      Input in(sq);
      in.with_weights(Weights{{FT(0), FT(0), FT(0), FT(0)}}).with_height(3);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 12.;
      e.top_area = 4.;
      e.rel_tol = 1e-6;
      run("prism (all vertical), H = 3", in, e);
    }

    {
      Input in(sq);
      in.with_weights(Weights{{FT(0), FT(0), FT(0), FT(0)}}).with_height(-3);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 12.;
      e.top_area = 4.;
      e.rel_tol = 1e-6;
      run("prism (all vertical), H = -3", in, e);
    }

    {
      // steep roof (weight 3): skeleton height 3, cut at 1.5 -> frustum of the square at offset 0.5
      Input in(sq);
      in.with_weights(Input::uniform(Polygon_with_holes_2(sq), 3.)).with_height(1.5);
      Expectations e;
      e.volume = 3 * (4 * 0.5 - 4 * 0.25 + 4 * 0.125 / 3);
      e.top_area = 1.;
      e.rel_tol = 1e-9;
      run("weight 3, H = 1.5", in, e);
    }

    {
      // angles as in the original test data
      Input in(SST::make_polygon<K>({{0,-5}, {0,0}, {5,0}, {5,5}, {-5,5}, {-5,-5}}));
      in.with_angles(Weights{{FT(90), FT(60), FT(60), FT(90), FT(60), FT(60)}}).with_height(6);
      Expectations e;
      e.volume = 162.37987499999997;
      e.rel_tol = 1e-5;
      run("polygon_000 / angles_000, H = 6", in, e);
    }
  }

  static void test_touching_offsets()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(1, 1, 3, 3, false));
    pwh.add_hole(SST::make_rectangle<K>(2, 6, 8, 8, false));
    Input in(pwh);
    in.with_height(1);
    Expectations e;
    e.genus = 2;
    run("offsets touching at the cut height", in, e);
  }

  static void run_all()
  {
    test_square_heights();
    test_rectangle_ridge_coincidence();
    test_ring_heights();
    test_multiple_holes_heights();
    test_height_with_weights();
    test_touching_offsets();
    test_event_height_sweeps();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  Max_height_tests<EPICK>::run_all();
  std::cout << "== EPECK ==" << std::endl;
  Max_height_tests<EPECK>::run_all();

  return SST::test_exit_code();
}
