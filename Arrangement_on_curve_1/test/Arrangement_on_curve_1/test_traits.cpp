// Tests the geometry traits that the other tests do not cover: the traits of
// lines in 2D, 3D, and dD, the adaptors of a 2D geometry traits with an open,
// a closed, a left-closed, a right-closed, and an identified parameter space,
// and the adaptor of a 2D geometry traits to the circle of directions.  Every permutation of
// every subset of a set of points is inserted, and then every vertex is
// removed on its own.  The order of the adaptor to the circle of directions is
// also compared with the order of `Circle_traits_1`.

#include <cstdlib>

#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/Line_2_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Line_3_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Line_d_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Open_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Closed_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Left_closed_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Right_closed_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Identified_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Direction_geom_traits_2_adaptor_1.h>

using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;

using Kernel_d = CGAL::Cartesian_d<CGAL::Quotient<CGAL::MP_Float>>;
using Point_d = Kernel_d::Point_d;

using Geom_traits_2 = CGAL::Arr_segment_traits_2<Kernel>;

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, std::size_t& failures) {
  std::size_t local_failures = 0;
  for_each_sequence(pts.size(), [&](const std::vector<int>& seq) {
    Arr arr(tr);
    std::vector<typename Arr::Point_1> inserted;
    for (int i : seq) {
      Aoc::insert(arr, pts[i]);
      inserted.push_back(pts[i]);
    }
    std::string why;
    if (! report(is_valid(arr, inserted, why), name + " (insert)", why, local_failures)) return;
    for (std::size_t k = 0; k < inserted.size(); ++k) {
      Arr copy(arr);
      copy.remove(find_vertex(copy, inserted[k]));
      auto rest = inserted;
      rest.erase(rest.begin() + k);
      report(is_valid(copy, rest, why), name + " (remove)", why, local_failures);
    }
  });
  failures += local_failures;
  std::cout << (local_failures ? "FAIL " : "PASS ") << name << std::endl;
}

// constructs the point of `Kernel_d` whose coordinates are `t` times (1, 2, 3).
Point_d point_d(int t) {
  std::vector<Kernel_d::FT> coords{Kernel_d::FT(t), Kernel_d::FT(2 * t), Kernel_d::FT(3 * t)};
  return Point_d(3, coords.begin(), coords.end());
}

int main() {
  std::size_t failures = 0;

  // Lines in 2D, 3D, and dD; the points are ordered along the direction of the line.
  {
    using Traits = Aoc::Line_2_traits_1<Kernel>;
    auto tr = std::make_shared<const Traits>(Kernel::Line_2(Point_2(0, 0), Point_2(1, 1)));
    test<Test_arrangement<Traits, false>>("Line_2_traits_1", tr,
                                          {Point_2(-3, -3), Point_2(1, 1), Point_2(2, 2), Point_2(5, 5)}, failures);
  }
  {
    using Traits = Aoc::Line_3_traits_1<Kernel>;
    auto tr = std::make_shared<const Traits>(Kernel::Line_3(Point_3(0, 0, 0), Point_3(1, 1, 1)));
    test<Test_arrangement<Traits, false>>("Line_3_traits_1", tr,
                                          {Point_3(-3, -3, -3), Point_3(1, 1, 1), Point_3(2, 2, 2), Point_3(5, 5, 5)},
                                          failures);
  }
  {
    using Traits = Aoc::Line_d_traits_1<Kernel_d>;
    std::vector<Kernel_d::FT> dir{Kernel_d::FT(1), Kernel_d::FT(2), Kernel_d::FT(3)};
    Kernel_d::Line_d line(point_d(0), Kernel_d::Direction_d(Kernel_d::Vector_d(3, dir.begin(), dir.end())));
    auto tr = std::make_shared<const Traits>(line);
    test<Test_arrangement<Traits, false>>("Line_d_traits_1", tr, {point_d(-1), point_d(1), point_d(2), point_d(4)},
                                          failures);
  }

  // Adaptors of a 2D geometry traits; the master curve is the segment from (0, 0) to (10, 5).
  auto tr_2 = std::make_shared<const Geom_traits_2>();
  const Point_2 left(0, 0), right(10, 5);
  {
    using Traits = Aoc::Open_geom_traits_2_adaptor_1<Geom_traits_2>;
    test<Test_arrangement<Traits, false>>("Open_geom_traits_2_adaptor_1", std::make_shared<const Traits>(tr_2),
                                          {Point_2(2, 1), Point_2(4, 2), Point_2(6, 3), Point_2(8, 4)}, failures);
  }
  {
    using Traits = Aoc::Closed_geom_traits_2_adaptor_1<Geom_traits_2>;
    test<Test_arrangement<Traits, false>>("Closed_geom_traits_2_adaptor_1",
                                          std::make_shared<const Traits>(tr_2, left, right),
                                          {left, Point_2(4, 2), Point_2(6, 3), right}, failures);
  }
  {
    using Traits = Aoc::Left_closed_geom_traits_2_adaptor_1<Geom_traits_2>;
    test<Test_arrangement<Traits, false>>("Left_closed_geom_traits_2_adaptor_1",
                                          std::make_shared<const Traits>(tr_2, left),
                                          {left, Point_2(4, 2), Point_2(6, 3), Point_2(8, 4)}, failures);
  }
  {
    using Traits = Aoc::Right_closed_geom_traits_2_adaptor_1<Geom_traits_2>;
    test<Test_arrangement<Traits, false>>("Right_closed_geom_traits_2_adaptor_1",
                                          std::make_shared<const Traits>(tr_2, right),
                                          {Point_2(2, 1), Point_2(4, 2), Point_2(6, 3), right}, failures);
  }

  // The adaptor of a 2D geometry traits with an identified parameter space;
  // the master curve is the equator of the sphere, which crosses the
  // identification curve of the sphere at (-1, 0, 0), the identification
  // point.  The points are ordered by their longitudes.
  {
    using Sphere_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
    using Traits = Aoc::Identified_geom_traits_2_adaptor_1<Sphere_traits_2>;
    auto sphere_tr = std::make_shared<const Sphere_traits_2>();
    auto ctr = sphere_tr->construct_point_2_object();
    auto tr = std::make_shared<const Traits>(sphere_tr);
    const std::vector<Traits::Point_1> pts{ctr(-1, 0, 0), ctr(1, 0, 0), ctr(0, 1, 0), ctr(0, -1, 0), ctr(1, 1, 0)};
    test<Test_arrangement<Traits, false>>("Identified_geom_traits_2_adaptor_1", tr, pts, failures);
    test<Test_arrangement<Traits, true>>("Identified_geom_traits_2_adaptor_1 (vector)", tr, pts, failures);

    // Only the point on the identification curve is the identification point,
    // and the other points are ordered by their longitudes.
    const auto is_on_id = tr->is_on_x_identification_1_object();
    const auto cmp = tr->compare_x_1_object();
    const bool ok = is_on_id(pts[0]) && ! is_on_id(pts[1]) && ! is_on_id(pts[2]) && ! is_on_id(pts[3]) &&
      (cmp(pts[3], pts[1]) == CGAL::SMALLER) && (cmp(pts[1], pts[4]) == CGAL::SMALLER) &&
      (cmp(pts[4], pts[2]) == CGAL::SMALLER) && (cmp(pts[2], pts[2]) == CGAL::EQUAL);
    report(ok, "Identified_geom_traits_2_adaptor_1 order", "order", failures);
    std::cout << (ok ? "PASS " : "FAIL ") << "Identified_geom_traits_2_adaptor_1 order" << std::endl;
  }

  // The adaptor to the circle of directions, of a kernel and of a 2D traits
  // derived from a kernel, with a reference direction on an axis and off the
  // axes; the directions include the reference direction and its opposite.
  {
    const std::vector<Direction_2> refs{Direction_2(1, 0), Direction_2(1, 2), Direction_2(-3, -1)};
    for (const auto& ref : refs) {
      std::vector<Direction_2> dirs{ref, -ref, Direction_2(1, 0), Direction_2(0, 1), Direction_2(-1, 0),
                                    Direction_2(0, -1), Direction_2(1, 1), Direction_2(2, 4), Direction_2(-1, 3)};
      const std::string suffix = " (reference " + std::to_string(&ref - refs.data()) + ")";
      {
        using Traits = Aoc::Direction_geom_traits_2_adaptor_1<Kernel>;
        test<Test_arrangement<Traits, false>>("Direction_geom_traits_2_adaptor_1<Kernel>" + suffix,
                                              std::make_shared<const Traits>(std::make_shared<const Kernel>(), ref),
                                              {ref, -ref, Direction_2(1, 1), Direction_2(-1, 3)}, failures);
        test<Test_arrangement<Traits, true>>("Direction_geom_traits_2_adaptor_1<Kernel> (vector)" + suffix,
                                             std::make_shared<const Traits>(Kernel(), ref),
                                             {ref, -ref, Direction_2(0, -1), Direction_2(2, -4)}, failures);
      }
      {
        using Traits = Aoc::Direction_geom_traits_2_adaptor_1<Geom_traits_2>;
        test<Test_arrangement<Traits, false>>("Direction_geom_traits_2_adaptor_1<Arr_segment_traits_2>" + suffix,
                                              std::make_shared<const Traits>(tr_2, ref),
                                              {ref, -ref, Direction_2(0, 1), Direction_2(1, 1)}, failures);
      }

      // The order, the identification point, and the equality agree with
      // those of `Circle_traits_1` with the same reference direction.
      using Traits = Aoc::Direction_geom_traits_2_adaptor_1<Kernel>;
      const Traits tr(Kernel(), ref);
      const Aoc::Circle_traits_1<Kernel> circle_tr(ref);
      const auto is_on_id = tr.is_on_x_identification_1_object();
      const auto circle_is_on_id = circle_tr.is_on_x_identification_1_object();
      const auto cmp = tr.compare_x_1_object();
      const auto circle_cmp = circle_tr.compare_x_1_object();
      const auto equal = tr.equal_1_object();
      bool ok = (tr.reference_direction() == ref);
      for (const auto& d1 : dirs) {
        ok = ok && (is_on_id(d1) == circle_is_on_id(d1));
        for (const auto& d2 : dirs) {
          ok = ok && (equal(d1, d2) == (d1 == d2));
          if (! is_on_id(d1) && ! is_on_id(d2)) ok = ok && (cmp(d1, d2) == circle_cmp(d1, d2));
        }
      }
      report(ok, "Direction_geom_traits_2_adaptor_1 vs. Circle_traits_1" + suffix, "order", failures);
      std::cout << (ok ? "PASS " : "FAIL ") << "Direction_geom_traits_2_adaptor_1 vs. Circle_traits_1" << suffix
                << std::endl;
    }
  }

  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
