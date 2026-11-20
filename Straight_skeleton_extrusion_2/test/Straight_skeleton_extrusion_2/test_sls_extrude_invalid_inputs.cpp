#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>

#include <iostream>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
struct Invalid_input_tests
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;
  using Outcome = SST::Outcome;

  // Check the outcome; when `false` is returned, the output must be empty
  static void expect(const std::string& name,
                     const Input& in,
                     const Outcome expected)
  {
    std::cout << "  test " << name << "..." << std::endl;

    const bool needs_checks = (expected == Outcome::PRECONDITION_EXCEPTION && !SST::has_cgal_preconditions) ||
                              (expected == Outcome::ASSERTION_EXCEPTION && !SST::has_cgal_assertions) ||
                              (expected == Outcome::POSTCONDITION_EXCEPTION && !SST::has_cgal_postconditions);
    if(needs_checks)
    {
      std::cout << "Skipping '" << name << "' (CGAL checks disabled)" << std::endl;
      return;
    }

    Mesh m;
    const Outcome o = extrude(in, m);
    CGAL_SLS_CHECK_MSG(o == expected, name << ": got " << to_string(o) << ", expected " << to_string(expected));
    if(o == Outcome::RETURNED_FALSE)
      CGAL_SLS_CHECK_MSG(num_vertices(m) == 0 && num_faces(m) == 0, name << ": output not empty after failure");
  }

  // Degenerate polygons
  static void test_degenerate_polygons()
  {
    std::cout << "test degenerate inputs" << std::endl;

    expect("empty polygon", Input(Polygon_2()), Outcome::RETURNED_FALSE);
    expect("single point", Input(SST::make_polygon<K>({{0,0}})), Outcome::RETURNED_FALSE);
    expect("two points", Input(SST::make_polygon<K>({{0,0}, {1,0}})), Outcome::RETURNED_FALSE);

    // Self-intersecting polygons are rejected by the skeleton construction
    expect("three collinear points", Input(SST::make_polygon<K>({{0,0}, {1,0}, {2,0}})), Outcome::RETURNED_FALSE);
    expect("bowtie", Input(SST::make_polygon<K>({{0,0}, {2,2}, {2,0}, {0,2}})), Outcome::RETURNED_FALSE);
    expect("self-overlapping", Input(SST::make_polygon<K>({{0,0}, {4,0}, {4,2}, {1,2}, {1,-1}, {0,-1}})), Outcome::RETURNED_FALSE);

    Input dup(SST::make_polygon<K>({{0,0}, {2,0}, {2,0}, {2,2}, {0,2}}));
    Mesh m;
    const Outcome o = extrude(dup, m);
    CGAL_SLS_CHECK_MSG(o == Outcome::SUCCESS || o == Outcome::RETURNED_FALSE, "duplicated vertex: " << to_string(o));

    // Weakly simple polygon with a zero-width antenna (repeated non-consecutive vertex)
    expect("antenna", Input(SST::make_polygon<K>({{0,0}, {4,0}, {4,4}, {2,4}, {2,6}, {2,4}, {0,4}})), Outcome::RETURNED_FALSE);
  }

  // Orientation
  static void test_orientation()
  {
    std::cout << "test orientation" << std::endl;

    Polygon_2 cw = SST::make_rectangle<K>(0, 0, 2, 2);
    cw.reverse_orientation();
    expect("clockwise outer boundary", Input(cw), Outcome::RETURNED_FALSE);

    Input in(cw);
    in.with_height(0.5);
    expect("clockwise outer boundary, H = 0.5", in, Outcome::RETURNED_FALSE);

    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
    pwh.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, true /*ccw: wrong*/));
    expect("counterclockwise hole", Input(pwh), Outcome::RETURNED_FALSE);
    Input out(pwh);
    out.with_weights(Input::uniform(pwh, -1.)).with_height(0.5);
    expect("counterclockwise hole, outward", out, Outcome::RETURNED_FALSE);
  }

  // Holes violating the "do not intersect" precondition: currently cleanly rejected
  static void test_touching_holes()
  {
    std::cout << "test touching" << std::endl;

    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
      pwh.add_hole(SST::make_rectangle<K>(0, 4, 2, 6, false)); // shares part of the edge x = 0
      expect("hole touching the outer boundary along an edge", Input(pwh), Outcome::RETURNED_FALSE);
    }

    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
      Polygon_2 h = SST::make_polygon<K>({{0,5}, {2,4}, {2,6}}); // vertex on the edge x = 0
      h.reverse_orientation();
      pwh.add_hole(h);
      expect("hole touching the outer boundary at a vertex", Input(pwh), Outcome::RETURNED_FALSE);
    }

    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
      pwh.add_hole(SST::make_rectangle<K>(2, 2, 4, 8, false));
      pwh.add_hole(SST::make_rectangle<K>(4, 2, 6, 8, false)); // shared edge x = 4
      expect("two holes sharing an edge", Input(pwh), Outcome::RETURNED_FALSE);
    }

    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 10));
      pwh.add_hole(SST::make_rectangle<K>(2, 2, 6, 6, false));
      pwh.add_hole(SST::make_rectangle<K>(4, 4, 8, 8, false)); // overlapping
      expect("overlapping holes", Input(pwh), Outcome::RETURNED_FALSE);
    }
  }

  // Invalid weights / angles / heights
  static void test_invalid_parameters()
  {
    std::cout << "test invalid params" << std::endl;

    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);

    {
      Input in(sq);
      in.with_weights(Weights{{FT(1), FT(-1), FT(1), FT(1)}});
      expect("mixed-sign weights", in, Outcome::RETURNED_FALSE);
      in.with_height(1);
      expect("mixed-sign weights, H = 1", in, Outcome::RETURNED_FALSE);
    }

    {
      Input in(sq);
      in.with_angles(Weights{{FT(45), FT(135), FT(45), FT(45)}});
      expect("mixed acute/obtuse angles", in, Outcome::RETURNED_FALSE);
    }

    {
      Input in(sq);
      in.with_weights(Weights{{FT(0), FT(0), FT(0), FT(0)}});
      expect("all vertical, no height", in, Outcome::RETURNED_FALSE);
    }

    {
      Input in(sq);
      in.with_weights(Weights{{FT(-1), FT(-1), FT(-1), FT(-1)}});
      expect("outward, no height", in, Outcome::RETURNED_FALSE);
    }

    // documented: "This parameter should not be null"
    {
      Input in(sq);
      in.with_height(0);
      expect("zero height", in, Outcome::PRECONDITION_EXCEPTION);
    }

    // documented: angles strictly within 0 and 180 degrees
    for(double a : {0., 180., -10., 200.})
    {
      Input in(sq);
      in.with_angles(Weights{{FT(a), FT(45), FT(45), FT(45)}});
      expect("angle " + std::to_string(a), in, Outcome::PRECONDITION_EXCEPTION);
    }

    // Collinear consecutive edges with different weights: the weighted skeleton construction itself
    // fails (create_interior_weighted_straight_skeleton_2() returns null), so `false` is returned.
    // Not documented for extrude_skeleton(). With equal weights, collinear edges are fine.
    {
      Input in(SST::make_polygon<K>({{0,0}, {1,0}, {2,0}, {2,2}, {0,2}}));
      in.with_weights(Weights{{FT(1), FT(1), FT(2), FT(1), FT(1)}});
      expect("collinear consecutive edges, different weights", in, Outcome::RETURNED_FALSE);
      in.with_weights(Weights{{FT(1), FT(2), FT(2), FT(1), FT(1)}});
      expect("collinear consecutive edges, equal weights", in, Outcome::SUCCESS);
    }

    // number of weights != number of edges
    {
      Input in(sq);
      in.with_weights(Weights{{FT(1), FT(1), FT(1)}});
      expect("too few weights", in, Outcome::PRECONDITION_EXCEPTION);
      in.with_weights(Weights{{FT(1), FT(1), FT(1), FT(1), FT(1)}});
      expect("too many weights", in, Outcome::PRECONDITION_EXCEPTION);
    }
  }

  static void test_vertex_point_map()
  {
    std::cout << "test out VPM" << std::endl;

    Mesh m;
    auto alt_vpm = m.template add_property_map<typename Mesh::Vertex_index, Point_3>("v:alt_point", Point_3(-42, -42, -42)).first;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(SST::make_rectangle<K>(0, 0, 2, 2),
                     m, CGAL::parameters::vertex_point_map(alt_vpm)));
    CGAL_SLS_CHECK(num_vertices(m) == 5);

    std::size_t filled = 0;
    for(auto v : vertices(m))
      if(get(alt_vpm, v) != Point_3(-42, -42, -42))
        ++filled;

    CGAL_SLS_CHECK_MSG(filled == num_vertices(m), "vpm ignored?");
  }

  static void test_sequence_container_input()
  {
    std::cout << "test sequence container" << std::endl;

    std::vector<Point_2> pts = {Point_2(0,0), Point_2(2,0), Point_2(2,2), Point_2(0,2)};
    Mesh m;
    CGAL_SLS_CHECK(CGAL::extrude_skeleton(pts, m));
    CGAL_SLS_CHECK(num_vertices(m) == 5);
  }

  static void run_all()
  {
    test_degenerate_polygons();
    test_orientation();
    test_touching_holes();
    test_invalid_parameters();
    test_vertex_point_map();
    test_sequence_container_input();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // the default, made explicit: failed preconditions / assertions throw
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  // CGAL_warning() is used on the clockwise + maximum height path: do not abort
  CGAL::set_warning_behaviour(CGAL::CONTINUE);

  std::cout << "== EPICK ==" << std::endl;
  Invalid_input_tests<EPICK>::run_all();
  std::cout << "== EPECK ==" << std::endl;
  Invalid_input_tests<EPECK>::run_all();

  return SST::test_exit_code();
}
