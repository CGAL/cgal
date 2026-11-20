#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
struct Outward_tests
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

  // For rectilinear inputs, the exterior offset at distance t of a convex rectilinear polygon has area:
  //   A + P*t + 4*t^2, (A is the base area, P is the perimeter)
  // and an a*b rectangular hole has area:
  //   (a - 2*t)(b - 2*t)
  // until it vanishes at:
  //   t = min(a, b) / 2.
  // With a slope of w, the offset distance at height z is z / |w|.

  // Integral over [0, H] of the area of the a x b rectangle shrinking at unit speed
  static double shrinking_rectangle_volume(double a, double b, double H)
  {
    const double T = (std::min)(H, 0.5 * (std::min)(a, b)); // see hole's max height above
    // (a - 2t)(b - 2t) = ab - 2(a+b) t + 4 t^2
    return a*b*T - (a + b)*T*T + 4.*T*T*T/3.;
  }

  // Integral over [0, H] of the area of the a x b rectangle growing at unit speed
  static double growing_rectangle_volume(double a, double b, double H)
  {
    return a*b*H + (a + b)*H*H + 4.*H*H*H/3.;
  }

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

  // Square [0,2]^2, weights -1: an inverted frustum, the top cap is the square grown by H
  static void test_square()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    for(double H : {1e-9, 0.5, 1., 2., 10., 1e6})
    {
      Input in(sq);
      in.with_weights(Input::uniform(in.pwh, -1.)).with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.genus = 0;
      e.zmax_abs = H;
      e.volume = growing_rectangle_volume(2, 2, H);
      e.rel_tol = 1e-9;
      Mesh m = run("outward square, H = " + std::to_string(H), in, e);

      if constexpr(SSE::internal::is_exact_kernel<K>)
      {
        if(H == 0.5) // exact: 4H + 4H^2 + 4H^3/3 = 2 + 1 + 1/6
          CGAL_SLS_CHECK(CGAL::Polygon_mesh_processing::volume(m) == FT(19) / FT(6));
      }
    }
  }

  // Obtuse angles are equivalent to negative weights: angle 135 == weight -1, 120 == -sqrt(3)
  static void test_angles()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    Input ina(sq);
    ina.with_angles(Input::uniform(ina.pwh, 135.)).with_height(1);
    Input inw(sq);
    inw.with_weights(Input::uniform(inw.pwh, -1.)).with_height(1);
    Mesh ma = run("outward square, angles 135", ina);
    Mesh mw = run("outward square, weights -1", inw);

    Input in120(sq);
    in120.with_angles(Input::uniform(in120.pwh, 120.)).with_height(std::sqrt(3.));
    Expectations e;
    e.volume = std::sqrt(3.) * growing_rectangle_volume(2, 2, 1.);
    e.top_area = 16.;
    e.rel_tol = 1e-9;
    run("outward square, angles 120", in120, e);
  }

  // Different weights per edge, and vertical walls mixed with outward slopes
  static void test_mixed_weights_and_vertical_walls()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    {
      // edges x = 0 and x = 2 (indices 0 and 2) go out at 45 degrees, y = 0 and y = 2 are vertical:
      // cross-section (2 + 2t) x 2, volume 6 at H = 1
      Input in(sq);
      in.with_angles(Weights{{FT(135), FT(90), FT(135), FT(90)}}).with_height(1);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 6.;
      e.top_area = 8.;
      e.rel_tol = 1e-6;
      run("outward with two vertical walls", in, e);
    }
    {
      // weights -1 on x edges, -2 on y edges: cross-section (2 + 2t)(2 + t), H = 1
      Input in(sq);
      in.with_weights(Weights{{FT(-1), FT(-2), FT(-1), FT(-2)}}).with_height(1);
      Expectations e;
      e.volume = 4. + 3. + 2. / 3.;
      e.top_area = 12.;
      e.rel_tol = 1e-9;
      run("outward with weights -1/-2", in, e);
    }
  }

  // Non-convex outer boundary: L-shape. Its exterior offset has area 12 + 16t + 4t^2
  // (5 convex corners add t^2 each, the reflex one removes t^2) and no event.
  static void test_L_shape()
  {
    const Polygon_2 L = SST::make_polygon<K>({{0,0}, {4,0}, {4,2}, {2,2}, {2,4}, {0,4}});
    for(double H : {1., 3.})
    {
      Input in(L);
      in.with_weights(Input::uniform(in.pwh, -1.)).with_height(H);
      Expectations e;
      e.genus = 0;
      e.volume = 12 * H + 8 * H * H + 4 * H * H * H / 3;
      e.rel_tol = 1e-9;
      e.top_area = 12 + 16 * H + 4 * H * H;
      run("outward L, H = " + std::to_string(H), in, e);
    }
  }

  // U-shape: the notch (width 2) is filled at t = 1 by the exterior offset (edge event of the exterior
  // skeleton). Cut before, at, and after the event.
  static void test_U_shape_events()
  {
    const Polygon_2 U = SST::make_polygon<K>({{0,0}, {6,0}, {6,4}, {4,4}, {4,2}, {2,2}, {2,4}, {0,4}});
    for(double H : {0.5, 1. - 1e-9, 1., 1. + 1e-9, 2.})
    {
      Input in(U);
      in.with_weights(Input::uniform(in.pwh, -1.)).with_height(H);
      Expectations e;
      e.genus = 0;
      e.check_self_intersections = (H != 1.);
      e.rel_tol = 1e-7;
      run("outward U, H = " + std::to_string(H), in, e);
    }
  }

  // Square ring [0,10]^2 \ [4,6]^2: the hole closes at t = 1.
  static void test_ring()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));

    for(double H : {0.5, 1. - 1e-9, 1., 1. + 1e-9, 2.})
    {
      Input in(pwh);
      in.with_weights(Input::uniform(pwh, -1.)).with_height(H);
      Expectations e;
      e.volume = growing_rectangle_volume(10, 10, H) - shrinking_rectangle_volume(2, 2, H);
      const double hole = (std::max)(0., 2 - 2 * H);
      e.top_area = (10 + 2 * H) * (10 + 2 * H) - hole * hole;
      e.rel_tol = 1e-7;
      if(H < 1) // at H == 1, the hole closes exactly at the top cap
        e.genus = 1; // the hole is still open at the top
      else if(H > 1)
        e.genus = 0; // the hole is closed by a pyramid-shaped 'ceiling'
      e.check_self_intersections = (H != 1.);
      run("outward ring, H = " + std::to_string(H), in, e);
    }
  }

  // Outer boundary growing, hole with vertical walls (a straight shaft through the solid)
  static void test_vertical_hole()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));
    Input in(pwh);
    in.with_weights(Weights{{FT(-1), FT(-1), FT(-1), FT(-1)}, {FT(0), FT(0), FT(0), FT(0)}}).with_height(2);
    Expectations e;
    e.genus = 1;
    e.volume = growing_rectangle_volume(10, 10, 2) - 4 * 2;
    e.rel_tol = 1e-6;
    e.top_area = 14 * 14 - 4;
    run("outward with vertical hole", in, e);
  }

  // Three holes closing at t = 1, 1.5 and 1
  static void test_multiple_holes()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(1, 1, 3, 3, false));
    pwh.add_hole(SST::make_rectangle<K>(6, 1, 9, 4, false));
    pwh.add_hole(SST::make_rectangle<K>(2, 6, 8, 8, false));

    for(double H : {0.25, 0.5, 0.999, 1.001, 1.25, 1.5, 2.})
    {
      Input in(pwh);
      in.with_weights(Input::uniform(pwh, -1.)).with_height(H);
      Expectations e;
      e.volume = growing_rectangle_volume(10, 10, H) - shrinking_rectangle_volume(2, 2, H)
               - shrinking_rectangle_volume(3, 3, H) - shrinking_rectangle_volume(6, 2, H);
      e.rel_tol = 1e-7;
      e.genus = (H < 1) ? 3 : (H < 1.5 ? 1 : 0);
      // at H = 1.5 the 3x3 hole closes exactly in the top cap (a single point): geometric
      // non-manifoldness, for which the documentation allows self-intersections
      e.check_self_intersections = (H != 1.5);
      run("outward three holes, H = " + std::to_string(H), in, e);
    }
  }

  // maximum_height(-H): mirror image (outward and downward)
  static void test_negative_heights()
  {
    Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0,0}, {6,0}, {6,4}, {4,4}, {4,2}, {2,2}, {2,4}, {0,4}}));
    Polygon_with_holes_2 ring(SST::make_rectangle<K>(0, 0, 10, 10));
    ring.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));

    for(const auto& p : {pwh, ring})
    {
      for(double H : {0.5, 2.})
      {
        Input up(p);
        up.with_weights(Input::uniform(p, -1.)).with_height(H);
        Input down(p);
        down.with_weights(Input::uniform(p, -1.)).with_height(-H);
        Mesh mu = run("outward up, H = " + std::to_string(H), up);
        Mesh md = run("outward down, H = " + std::to_string(-H), down);
        SST::flip_mesh(md);
        CGAL_SLS_CHECK(PMP::are_equal(mu, md));
      }
    }
  }

  // Generic non-convex input with non-axis-aligned edges and random-ish weights
  static void test_generic()
  {
    Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0,3}, {-1,1}, {-3,1}, {-1.5,-0.5}, {-2,-3},
                                             {0,-1.5}, {2,-3}, {1.5,-0.5}, {3,1}, {1,1}}));
    Polygon_2 hole = SST::make_polygon<K>({{-0.5,-0.25}, {0.5,-0.25}, {0,0.5}});
    hole.reverse_orientation();
    pwh.add_hole(hole);

    for(double H : {0.1, 0.5, 2.})
    {
      Input in(pwh);
      // non-uniform weights on the outer boundary; uniform on the hole
      in.with_weights(Weights{{FT(-1), FT(-2), FT(-1.5), FT(-1), FT(-3), FT(-1), FT(-1), FT(-2), FT(-1), FT(-1)},
                              {FT(-2), FT(-2), FT(-2)}}).with_height(H);
      Expectations e;
      run("outward star with hole, H = " + std::to_string(H), in, e);
    }
  }

  // An outward extrusion requires a specified maximum height.
  static void test_requires_height()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    Input in(sq);
    in.with_weights(Input::uniform(in.pwh, -1.));
    Mesh m;
    CGAL_SLS_CHECK(extrude(in, m) == Outcome::RETURNED_FALSE);
    CGAL_SLS_CHECK(num_faces(m) == 0);

    in.with_height(0);
    CGAL_SLS_CHECK(extrude(in, m) == (SST::has_cgal_preconditions ? Outcome::PRECONDITION_EXCEPTION
                                                                   : Outcome::RETURNED_FALSE));

    Input ina(sq);
    ina.with_angles(Input::uniform(ina.pwh, 100.));
    CGAL_SLS_CHECK(extrude(ina, m) == Outcome::RETURNED_FALSE);
  }

  // to check that holes have the correct weights
  static void test_outward_hole_weights()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 100, 10));
    pwh.add_hole(SST::make_rectangle<K>(3, 3, 70, 7, false));

    Input inw(pwh);
    inw.with_weights(Weights{{FT(1), FT(1), FT(1), FT(1)}, {FT(1), FT(10), FT(100), FT(1)}});
    run("inward, non-uniform hole weights", inw);

    Input out(pwh);
    out.with_weights(Weights{{FT(-1), FT(-1), FT(-1), FT(-1)}, {FT(-1), FT(-10), FT(-100), FT(-1)}}).with_height(0.5);
    run("outward, non-uniform hole weights", out);
  }

  static void run_all()
  {
    test_outward_hole_weights();
    test_square();
    test_angles();
    test_mixed_weights_and_vertical_walls();
    test_L_shape();
    test_U_shape_events();
    test_ring();
    test_vertical_hole();
    test_multiple_holes();
    test_negative_heights();
    test_generic();
    test_requires_height();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  Outward_tests<EPICK>::run_all();
  std::cout << "== EPECK ==" << std::endl;
  Outward_tests<EPECK>::run_all();

  return SST::test_exit_code();
}
