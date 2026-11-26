#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_batched_point_location.h>
#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Default_visitor.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/Double_grid_snap_rounding_traits_2.h>
#include <CGAL/Named_function_parameters.h>

#include <atomic>
#include <vector>
#include <cassert>
#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Point_2 = Traits::Point_2;
using Segment_2 = Traits::X_monotone_curve_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits>;
using Overlay_traits = CGAL::Arr_default_overlay_traits<Arrangement_2>;

template <typename GeometryTraits_2>
class Test_sweep_visitor : public CGAL::Surface_sweep_2::Default_visitor<
                             Test_sweep_visitor<GeometryTraits_2>,
                             GeometryTraits_2> {
public:
  std::size_t event_count{0};

  template <typename Event>
  bool before_handle_event(Event*) {
    ++event_count;
    return true;
  }
};

void test_named_parameters() {
  std::cout << "Testing named parameter stop_predicate..." << std::endl;
  std::atomic<bool> flag{false};

  auto np = CGAL::parameters::stop_predicate([&] { return flag.load(); });
  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;

  auto pred = choose_parameter(get_parameter(np, CGAL::internal_np::stop_predicate), std::function<bool()>());
  std::function<bool()> fn = pred;
  assert(bool(fn));
  assert(!fn());

  flag.store(true);
  assert(fn());

  // Default value should be empty std::function
  auto np_default = CGAL::parameters::default_values();
  auto pred_default = choose_parameter(get_parameter(np_default, CGAL::internal_np::stop_predicate), std::function<bool()>());
  assert(!pred_default);
  std::cout << "  -> Named parameters test PASSED" << std::endl;
}

void test_surface_sweep_cancellation() {
  std::cout << "Testing Surface_sweep_2 cancellation..." << std::endl;

  Traits traits;
  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0, 0), Point_2(4, 0)),
    Segment_2(Point_2(0, 1), Point_2(4, 1)),
    Segment_2(Point_2(0, 2), Point_2(4, 2)),
    Segment_2(Point_2(0, 3), Point_2(4, 3))
  };

  // Case 1: Normal sweep completion
  {
    Test_sweep_visitor<Traits> visitor;
    CGAL::Surface_sweep_2::No_intersection_surface_sweep_2<Test_sweep_visitor<Traits>> sweep(&traits, &visitor);
    sweep.sweep(segments.begin(), segments.end());
    assert(sweep.is_completed());
    assert(!sweep.is_stopped());
    assert(visitor.event_count > 0);
  }

  // Case 2: Early cancellation before sweep starts
  {
    Test_sweep_visitor<Traits> visitor;
    CGAL::Surface_sweep_2::No_intersection_surface_sweep_2<Test_sweep_visitor<Traits>> sweep(&traits, &visitor);
    sweep.set_stop_predicate([] { return true; });
    sweep.sweep(segments.begin(), segments.end());
    assert(sweep.is_stopped());
    assert(!sweep.is_completed());
    assert(visitor.event_count == 0);
  }

  // Case 3: Mid-sweep cancellation after 2 events
  {
    Test_sweep_visitor<Traits> visitor;
    CGAL::Surface_sweep_2::No_intersection_surface_sweep_2<Test_sweep_visitor<Traits>> sweep(&traits, &visitor);
    std::size_t checks = 0;
    sweep.set_stop_predicate([&checks] {
      return (++checks >= 2);
    });
    sweep.sweep(segments.begin(), segments.end());
    assert(sweep.is_stopped());
    assert(!sweep.is_completed());
  }

  std::cout << "  -> Surface_sweep_2 cancellation test PASSED" << std::endl;
}

void test_overlay_cancellation() {
  std::cout << "Testing CGAL::overlay cancellation..." << std::endl;

  Arrangement_2 arr1, arr2;

  // Box 1: [0, 2] x [0, 2]
  std::vector<Segment_2> segs1 = {
    Segment_2(Point_2(0, 0), Point_2(2, 0)),
    Segment_2(Point_2(2, 0), Point_2(2, 2)),
    Segment_2(Point_2(0, 2), Point_2(2, 2)),
    Segment_2(Point_2(0, 0), Point_2(0, 2))
  };
  CGAL::insert_non_intersecting_curves(arr1, segs1.begin(), segs1.end());

  // Box 2: [1, 3] x [1, 3]
  std::vector<Segment_2> segs2 = {
    Segment_2(Point_2(1, 1), Point_2(3, 1)),
    Segment_2(Point_2(3, 1), Point_2(3, 3)),
    Segment_2(Point_2(1, 3), Point_2(3, 3)),
    Segment_2(Point_2(1, 1), Point_2(1, 3))
  };
  CGAL::insert_non_intersecting_curves(arr2, segs2.begin(), segs2.end());

  Overlay_traits ot;

  // Case 1: Overlay completes successfully
  {
    Arrangement_2 res;
    bool completed = CGAL::overlay(
      arr1, arr2, res, ot,
      CGAL::parameters::stop_predicate([] { return false; }));
    assert(completed);
    assert(!res.is_empty());
    assert(res.number_of_vertices() > 0);
  }

  // Case 2: Overlay is cancelled via stop_predicate
  {
    Arrangement_2 res;
    std::atomic<bool> cancelled{true};
    bool completed = CGAL::overlay(
      arr1, arr2, res, ot,
      CGAL::parameters::stop_predicate([&] { return cancelled.load(); }));
    assert(!completed);
    assert(res.is_empty());
  }

  std::cout << "  -> CGAL::overlay cancellation test PASSED" << std::endl;
}

void test_insertion_cancellation() {
  std::cout << "Testing CGAL::insert and insert_non_intersecting_curves cancellation..." << std::endl;

  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0, 0), Point_2(4, 0)),
    Segment_2(Point_2(0, 1), Point_2(4, 1)),
    Segment_2(Point_2(0, 2), Point_2(4, 2))
  };

  // Case 1: Successful non-intersecting insertion
  {
    Arrangement_2 arr;
    bool completed = CGAL::insert_non_intersecting_curves(
      arr, segments.begin(), segments.end(),
      CGAL::parameters::stop_predicate([] { return false; }));
    assert(completed);
    assert(!arr.is_empty());
    assert(arr.number_of_edges() == 3);
  }

  // Case 2: Cancelled non-intersecting insertion
  {
    Arrangement_2 arr;
    bool completed = CGAL::insert_non_intersecting_curves(
      arr, segments.begin(), segments.end(),
      CGAL::parameters::stop_predicate([] { return true; }));
    assert(!completed);
    assert(arr.is_empty());
  }

  // Intersecting segments:
  std::vector<Segment_2> x_segs = {
    Segment_2(Point_2(0, 0), Point_2(4, 4)),
    Segment_2(Point_2(0, 4), Point_2(4, 0))
  };

  // Case 3: Successful general insertion
  {
    Arrangement_2 arr;
    bool completed = CGAL::insert(
      arr, x_segs.begin(), x_segs.end(),
      CGAL::parameters::stop_predicate([] { return false; }));
    assert(completed);
    assert(!arr.is_empty());
  }

  // Case 4: Cancelled general insertion
  {
    Arrangement_2 arr;
    bool completed = CGAL::insert(
      arr, x_segs.begin(), x_segs.end(),
      CGAL::parameters::stop_predicate([] { return true; }));
    assert(!completed);
    assert(arr.is_empty());
  }

  std::cout << "  -> CGAL::insert cancellation test PASSED" << std::endl;
}

void test_locate_with_stop_predicate() {
  std::cout << "Testing CGAL::locate with stop_predicate..." << std::endl;

  Arrangement_2 arr;
  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0, 0), Point_2(4, 0)),
    Segment_2(Point_2(0, 2), Point_2(4, 2))
  };
  CGAL::insert_non_intersecting_curves(arr, segments.begin(), segments.end());

  std::vector<Point_2> query_points = { Point_2(1, 1), Point_2(2, 0) };
  using Point_location_result = CGAL::Arr_point_location_result<Arrangement_2>;
  using Query_result = std::pair<Point_2, Point_location_result::Type>;
  std::vector<Query_result> results;

  CGAL::locate(arr, query_points.begin(), query_points.end(),
               std::back_inserter(results),
               CGAL::parameters::stop_predicate([] { return false; }));
  assert(results.size() == 2);

  std::cout << "  -> CGAL::locate test PASSED" << std::endl;
}

void test_snap_rounding_cancellation() {
  std::cout << "Testing CGAL::snap_rounding_2 cancellation..." << std::endl;

  using Seg = Kernel::Segment_2;
  using Pt = Kernel::Point_2;
  using Polyline_2 = std::vector<Pt>;

  std::vector<Seg> segs = {
    Seg(Pt(0, 0), Pt(2, 2)),
    Seg(Pt(0, 2), Pt(2, 0))
  };

  // Case 1: Normal snap rounding
  {
    std::vector<Polyline_2> out;
    CGAL::snap_rounding_2(
      segs, out,
      CGAL::parameters::geom_traits(CGAL::Double_grid_snap_rounding_traits_2<Kernel>())
                       .stop_predicate([] { return false; }));
    assert(!out.empty());
  }

  // Case 2: Cancelled snap rounding
  {
    std::vector<Polyline_2> out;
    CGAL::snap_rounding_2(
      segs, out,
      CGAL::parameters::geom_traits(CGAL::Double_grid_snap_rounding_traits_2<Kernel>())
                       .stop_predicate([] { return true; }));
    assert(out.empty());
  }

  std::cout << "  -> CGAL::snap_rounding_2 cancellation test PASSED" << std::endl;
}

int main() {
  std::cout << "========================================" << std::endl;
  std::cout << "Running CGAL stop_predicate Test Suite" << std::endl;
  std::cout << "========================================" << std::endl;

  test_named_parameters();
  test_surface_sweep_cancellation();
  test_overlay_cancellation();
  test_insertion_cancellation();
  test_locate_with_stop_predicate();
  test_snap_rounding_cancellation();

  std::cout << "========================================" << std::endl;
  std::cout << "ALL 6 TEST SUITES PASSED SUCCESSFULLY!" << std::endl;
  std::cout << "========================================" << std::endl;

  return 0;
}
