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
struct Inward_tests
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

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

  // Square [0,2]^2, uniform weights: a pyramid with apex (1,1,1).
  // 5 vertices; 4 roof triangles + 2 bottom triangles; volume = 4 / 3.
  static void test_square_pyramid()
  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));

    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.genus = 0;
    e.volume = 4. / 3.;
    e.zmax_abs = 1.;

    Mesh m = run("square pyramid", in, e);

    if constexpr(SSE::internal::is_exact_kernel<K>)
    {
      CGAL_SLS_CHECK(SST::has_vertex(m, Point_3(1, 1, 1)));
      CGAL_SLS_CHECK(CGAL::Polygon_mesh_processing::volume(m) == FT(4) / FT(3));
    }
    else
    {
      CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(1, 1, 1)));
    }
  }

  // Rectangle [0,4]x[0,2]: hip roof, ridge from (1,1,1) to (3,1,1).
  // 6 vertices; 2 triangular + 2 trapezoidal roof faces (2 triangles each) + 2 bottom = 8 faces.
  // Volume of a hip roof of width W, length L, slope 1: W^2*L / 4 - W^3 / 12 = 10/3.
  static void test_rectangle_hip_roof()
  {
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));

    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.genus = 0;
    e.volume = 10. / 3.;
    e.zmax_abs = 1.;

    Mesh m = run("rectangle hip roof", in, e);
    if constexpr(SSE::internal::is_exact_kernel<K>)
    {
      CGAL_SLS_CHECK(SST::has_vertex(m, Point_3(1, 1, 1)));
      CGAL_SLS_CHECK(SST::has_vertex(m, Point_3(3, 1, 1)));
      CGAL_SLS_CHECK(CGAL::Polygon_mesh_processing::volume(m) == FT(10) / FT(3));
    }
    else
    {
      CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(1, 1, 1)));
      CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(3, 1, 1)));
    }
  }

  // L-shape (one reflex vertex, legs of width 2): all skeleton events at z = 1.
  static void test_L_shape()
  {
    Input in(SST::make_polygon<K>({{0,0}, {4,0}, {4,2}, {2,2}, {2,4}, {0,4}}));

    Expectations e;
    e.nv = 9;
    e.nf = 14;
    e.genus = 0;
    e.volume = 16. / 3.;
    e.zmax_abs = 1.;

    run("L-shape", in, e);
  }

  // T-shape (two reflex vertices)
  static void test_T_shape()
  {
    Input in(SST::make_polygon<K>({{0,0}, {6,0}, {6,2}, {4,2}, {4,6}, {2,6}, {2,2}, {0,2}}));

    Expectations e;
    e.nv = 12;
    e.nf = 20;
    e.genus = 0;
    e.volume = 28. / 3.;
    e.zmax_abs = 1.;

    run("T-shape", in, e);
  }

  // Non-convex polygon with non axis-aligned edges and several reflex vertices (split events)
  static void test_star()
  {
    Input in(SST::make_polygon<K>({{0,3}, {-1,1}, {-3,1}, {-1.5,-0.5}, {-2,-3}, {0,-1.5}, {2,-3}, {1.5,-0.5}, {3,1}, {1,1}}));

    Expectations e;
    e.genus = 0;

    run("star", in, e);
  }

  // Uniform weights but != 1
  static void test_uniform_weight_scaling()
  {
    Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    Input in(sq);
    in.with_weights(Input::uniform(in.pwh, 2.));
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 8. / 3.;
    e.zmax_abs = 2.;
    run("square, weight 2", in, e);

    Input in_small(sq);
    in_small.with_weights(Input::uniform(in_small.pwh, 0.25));
    Expectations e2;
    e2.nv = 5;
    e2.nf = 6;
    e2.volume = 1. / 3.;
    e2.zmax_abs = 0.25;
    run("square, weight 1/4", in_small, e2);
  }

  // angles: 45 degrees == weight 1, 60 degrees == sqrt(3)
  static void test_angles()
  {
    Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);

    Input in45(sq);
    in45.with_angles(Input::uniform(in45.pwh, 45.));
    Expectations e45;
    e45.nv = 5;
    e45.nf = 6;
    e45.volume = 4. / 3.;
    e45.zmax_abs = 1.;
    e45.rel_tol = 1e-12;
    Mesh m45 = run("square, angles 45", in45, e45);
    Mesh m1 = run("square, weights 1", Input(sq), e45);

    const double tol = 1e-9 * SST::diag_xy(PMP::bbox(m45));
    CGAL_SLS_CHECK(SST::approx_equal(PMP::volume(m45), PMP::volume(m1), tol));

    Input in60(sq);
    in60.with_angles(Input::uniform(in60.pwh, 60.));
    Expectations e60;
    e60.nv = 5;
    e60.volume = 4. * std::sqrt(3.) / 3.;
    e60.zmax_abs = std::sqrt(3.);
    run("square, angles 60", in60, e60);
  }

  // `angles` should be ignored when `weights` is provided
  static void test_weights_over_angles()
  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_weights(Input::uniform(in.pwh, 2.));
    Weights angles = Input::uniform(in.pwh, 30.);

    Mesh m;
    bool res = CGAL::extrude_skeleton(in.pwh, m, CGAL::parameters::weights(in.weights).angles(angles));
    CGAL_SLS_CHECK(res);
    Expectations e;
    e.zmax_abs = 2.;
    e.volume = 8. / 3.;
    check_extrusion("weights + angles", in, m, e);
  }

  // The i-th weight is associated with edge (v_{i-1}, v_i)
  // Weights {1,2,1,1} on the square make the edge v0v1 (y = 0) twice as steep => the ridge
  // goes from (1, 1/2, 1) to (1, 1, 1).
  static void test_weight_edge_association()
  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_weights(Weights{{FT(1), FT(2), FT(1), FT(1)}});
    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.volume = 1.5;
    e.zmax_abs = 1.;
    Mesh m = run("weight edge association", in, e);
    CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(1, 0.5, 1)));
    CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(1, 1, 1)));
    if constexpr(SSE::internal::is_exact_kernel<K>)
      CGAL_SLS_CHECK(SST::has_vertex(m, Point_3(1, FT(1)/FT(2), 1)));
    else
      CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(1, FT(1)/FT(2), 1)));
  }

  // Different weights producing a single simultaneous (4-way) event: rectangle [0,4]x[0,2] whose short
  // sides have weight 1/2. Both pairs of opposite edges meet at z = 1: a pyramid with apex (2,1,1).
  static void test_weights_simultaneous_event()
  {
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));
    in.with_weights(Weights{{FT(0.5), FT(1), FT(0.5), FT(1)}});
    Expectations e;
    e.nv = 5;
    e.nf = 6;
    e.volume = 8. / 3.;
    e.zmax_abs = 1.;
    Mesh m = run("weighted simultaneous event", in, e);
    if constexpr(SSE::internal::is_exact_kernel<K>)
      CGAL_SLS_CHECK(SST::has_vertex(m, Point_3(2, 1, 1)));
    else
      CGAL_SLS_CHECK(SST::has_vertex_approx(m, Point_3(2, 1, 1)));
  }

  // vertical walls (angle=90 / weight=0)
  static void test_vertical_walls()
  {
    // one vertical edge (y = 0): region at height t is [t, 2-t] x [0, 2-t]; volume = 5/3
    {
      Input in(SST::make_rectangle<K>(0, 0, 2, 2));
      in.with_weights(Weights{{FT(1), FT(0), FT(1), FT(1)}});
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = 5. / 3.;
      e.zmax_abs = 1.;
      e.rel_tol = 1e-6;
      run("one vertical wall (weight 0)", in, e);
    }

    // gable roof on [0,4]x[0,2], ridge (2,0,2)-(2,2,2),
    // volume = triangle (base 4, height 2) * length 2 = 8.
    {
      Input in(SST::make_rectangle<K>(0, 0, 4, 2));
      in.with_angles(Weights{{FT(45), FT(90), FT(45), FT(90)}});
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = 8.;
      e.zmax_abs = 2.;
      e.rel_tol = 1e-6;
      run("gable roof (angles 45/90)", in, e);
    }

    // two adjacent vertical walls: region at height t is [t, 2] x [0, 2-t], so volume = 8/3
    {
      Input in(SST::make_rectangle<K>(0, 0, 2, 2));
      in.with_angles(Weights{{FT(45), FT(90), FT(90), FT(45)}});
      Expectations e;
      e.nv = 5;
      e.nf = 6;
      e.volume = 8. / 3.;
      e.zmax_abs = 2.;
      e.rel_tol = 1e-6;
      run("two adjacent vertical walls", in, e);
    }

    // three vertical walls: shed roof, region at height t is [0, 2] x [0, 2-t], so volume = 4
    {
      Input in(SST::make_rectangle<K>(0, 0, 2, 2));
      in.with_angles(Weights{{FT(90), FT(90), FT(90), FT(45)}});
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = 4.;
      e.zmax_abs = 2.;
      e.rel_tol = 1e-6;
      run("three vertical walls", in, e);
    }
  }

  // Square ring [0,10]^2 minus [4,6]^2: the inner and outer fronts meet at t = 2.
  // Cross-section area at height t: (10-2t)^2 - (2+2t)^2 = 96 - 48t, volume = 96. Genus 1.
  static void test_square_ring()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false /*cw*/));
    Input in(pwh);
    Expectations e;
    e.nv = 12;
    e.nf = 24;
    e.genus = 1;
    e.volume = 96.;
    e.zmax_abs = 2.;
    Mesh m = run("square ring", in, e);
    if constexpr(SSE::internal::is_exact_kernel<K>)
      CGAL_SLS_CHECK(CGAL::Polygon_mesh_processing::volume(m) == FT(96));
  }

  // Several holes of different sizes (holes consumed at different heights); per-contour weights
  static void test_multiple_holes()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(1, 1, 3, 3, false));
    pwh.add_hole(SST::make_rectangle<K>(6, 1, 9, 4, false));
    pwh.add_hole(SST::make_rectangle<K>(2, 6, 8, 8, false));

    Input in(pwh);
    Expectations e;
    e.genus = 3;
    e.zmax_abs = 1.5;
    run("three holes, uniform", in, e);

    Input inw(pwh);
    inw.with_weights(Weights{{FT(1), FT(2), FT(1), FT(2)},
                             {FT(3), FT(3), FT(3), FT(3)},
                             {FT(0.5), FT(0.5), FT(0.5), FT(0.5)},
                             {FT(1), FT(1.5), FT(2), FT(2.5)}});
    Expectations ew;
    ew.genus = 3;
    run("three holes, per-contour weights", inw, ew);
  }

  static void test_nonconvex_hole()
  {
    // non-convex hole inside a non-convex outer boundary
    Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0,0}, {12,0}, {12,8}, {7,8}, {7,12}, {0,12}}));
    Polygon_2 hole = SST::make_polygon<K>({{2,2}, {6,2}, {6,4}, {4,4}, {4,6}, {2,6}});
    hole.reverse_orientation();
    pwh.add_hole(hole);
    Input in(pwh);
    Expectations e;
    e.genus = 1;
    e.zmax_abs = 3;
    run("non-convex hole in non-convex polygon", in, e);
  }

  static void test_negative_heights()
  {
    Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0,0}, {6,0}, {6,2}, {4,2}, {4,6}, {0,6}}));

    Polygon_with_holes_2 ring(SST::make_rectangle<K>(0, 0, 10, 10));
    ring.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));

    for(const auto& p : {pwh, ring})
    {
      for(double H : {0.5, 1., 1.5, 5.})
      {
        Input up(p);
        up.with_height(H);
        Input down(p);
        down.with_height(-H);
        Expectations e;
        e.genus = int(p.number_of_holes());
        Mesh mu = run("upward H = " + std::to_string(H), up, e);
        Mesh md = run("downward H = " + std::to_string(-H), down, e);
        SST::flip_mesh(md);
        CGAL_SLS_CHECK_MSG(PMP::are_equal(mu, md), "mirror, H = " << H);
      }
    }

    // a huge negative height is not the "no height" sentinel, but must behave like no height
    Input down(pwh);
    down.with_height(-1e300);
    Expectations e;
    e.zmax_abs = 2.;
    run("downward, H = -1e300", down, e);
  }

  static void test_P_vs_PWH()
  {
    std::cout << "test P vs PWH" << std::endl;

    Polygon_2 p = SST::make_polygon<K>({{0,0}, {5,0}, {5,1}, {2,3}, {0,1}});
    Mesh m1, m2;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m1));
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(Polygon_with_holes_2(p), m2));
    CGAL_SLS_CHECK(PMP::are_equal(m1, m2));
  }

  static void test_polyhedron_output()
  {
    std::cout << "test polyhedron output" << std::endl;

    using Polyhedron = CGAL::Polyhedron_3<K>;
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));
    Polyhedron P;
    const Outcome o = extrude(in, P);
    CGAL_SLS_CHECK(o == Outcome::SUCCESS);
    Expectations e;
    e.nv = 6;
    e.nf = 8;
    e.volume = 10. / 3.;
    check_extrusion("Polyhedron_3 output", in, P, e);
  }

  static void test_NP()
  {
    std::cout << "test NP" << std::endl;

    Polygon_2 p = SST::make_rectangle<K>(0, 0, 4, 2);
    Weights w{{FT(1), FT(1), FT(1), FT(1)}};
    Mesh m1, m2;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m1, CGAL::parameters::weights(w).geom_traits(K())));
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m2));
    CGAL_SLS_CHECK(PMP::are_equal(m1, m2));

    // a const weight range is accepted too
    const Weights cw{{FT(2), FT(2), FT(2), FT(2)}};
    Mesh m3;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m3, CGAL::parameters::weights(cw)));
    CGAL_SLS_CHECK(SST::approx_equal(PMP::volume(m3), 2. * 10. / 3., 1e-12));
  }

  // Undocumented: the output mesh is not cleared, the new extrusion is appended as a new
  // connected component. This test pins down the current behavior.
  static void test_output_is_appended()
  {
    std::cout << "test output_is_appended" << std::endl;

    Polygon_2 p = SST::make_rectangle<K>(0, 0, 2, 2);
    Mesh m;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m));
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(p, m));
    CGAL_SLS_CHECK(num_vertices(m) == 10 && num_faces(m) == 12);
  }

  static void run_all()
  {
    test_square_pyramid();
    test_rectangle_hip_roof();
    test_L_shape();
    test_T_shape();
    test_star();
    test_uniform_weight_scaling();
    test_angles();
    test_weights_over_angles();
    test_weight_edge_association();
    test_weights_simultaneous_event();
    test_vertical_walls();
    test_square_ring();
    test_multiple_holes();
    test_nonconvex_hole();
    test_negative_heights();
    test_P_vs_PWH();
    test_polyhedron_output();
    test_NP();
    test_output_is_appended();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  Inward_tests<EPICK>::run_all();
  std::cout << "== EPECK ==" << std::endl;
  Inward_tests<EPECK>::run_all();

  return SST::test_exit_code();
}
