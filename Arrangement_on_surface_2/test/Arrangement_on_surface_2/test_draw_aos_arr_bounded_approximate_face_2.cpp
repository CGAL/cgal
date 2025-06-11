#include "CGAL/Bbox_2.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include <CGAL/Draw_aos/Arr_bounded_approximate_face_2.h>
#include <cstddef>
#include <iterator>
#include <vector>

using Approx_point = CGAL::Arr_approximation_geometry_traits::Approx_point;

struct Test_case
{
  std::vector<Approx_point> points_to_insert;
  std::vector<Approx_point> expected_points;
  CGAL::Bbox_2 bbox;

  Test_case(std::vector<Approx_point> points, std::vector<Approx_point> expected, CGAL::Bbox_2 b)
      : points_to_insert(std::move(points))
      , expected_points(std::move(expected))
      , bbox(b) {}
};

void test_geom_fill_corners_output_iterator(const Test_case& test_case) {
  std::vector<Approx_point> points{};
  auto inserter = std::back_inserter(points);
  CGAL::Geom_fill_corners_output_iterator<std::back_insert_iterator<decltype(points)>> geom_inserter(
      inserter, CGAL::Bbox_2(test_case.bbox));

  for(const auto& pt : test_case.points_to_insert) {
    geom_inserter++ = pt;
  }

  if(points != test_case.expected_points) {
    std::cerr << "Test failed for bbox: " << test_case.bbox << std::endl;
    std::cerr << "Points inserted: \n";
    for(const auto& pt : test_case.points_to_insert) {
      std::cerr << pt << "\n";
    }
    std::cerr << "\nExpected: \n";
    for(const auto& pt : test_case.expected_points) {
      std::cerr << pt << "\n";
    }
    std::cerr << "\nGot:\n";
    for(const auto& pt : points) {
      std::cerr << pt << "\n";
    }
    std::cerr << std::endl;
    assert(false);
  }
}

int main() {
  CGAL::Bbox_2 bbox(0, 0, 10, 10);
  std::vector<Test_case> cases{
      /**
       *  --------------------
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1-------------------
       */
      Test_case({Approx_point(0, 0)}, {Approx_point(0, 0)}, bbox),

      /**
       *  2-------------------
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1-------------------
       *
       *  Expected:
       *
       *  4------------------3
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1------------------2
       */
      Test_case({Approx_point(0, 0), Approx_point(0, 10)},
                {Approx_point(0, 0), Approx_point(10, 0), Approx_point(10, 10), Approx_point(0, 10)}, bbox),
      /**
       *  --------------------
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1------------------2
       *
       *  Expected:
       *
       *  2-------------------
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1------------------2
       */
      Test_case({Approx_point(0, 0), Approx_point(10, 0)}, {Approx_point(0, 0), Approx_point(10, 0)}, bbox),
      /**
       *  --------------------
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1---------2---------
       *
       *  Expected:
       *
       *  --------------------
       *  |                  |
       *  |                  |
       *  2                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  1-------------------
       */
      Test_case({Approx_point(0, 0), Approx_point(5, 0)}, {Approx_point(0, 0), Approx_point(5, 0)}, bbox),
      /**
       *  --------------------
       *  |                  |
       *  |                  |
       *  |                  1
       *  |                  |
       *  |                  2
       *  |                  |
       *  --------------------
       *
       *  Expected:
       *
       *  3------------------2
       *  |                  |
       *  |                  |
       *  |                  1
       *  |                  |
       *  |                  6
       *  |                  |
       *  4------------------5
       */
      Test_case({Approx_point(10, 5), Approx_point(10, 3)},
                {Approx_point(10, 5), Approx_point(10, 10), Approx_point(0, 10), Approx_point(0, 0),
                 Approx_point(10, 0), Approx_point(10, 3)},
                bbox),
      /**
       *  -------------------1
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  2
       *  |                  |
       *  -------------------|
       *
       *  Expected:
       *
       *  2------------------1
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  |
       *  |                  5
       *  |                  |
       *  3------------------4
       */
      Test_case({Approx_point(10, 10), Approx_point(10, 3)},
                {
                    Approx_point(10, 10),
                    Approx_point(0, 10),
                    Approx_point(0, 0),
                    Approx_point(10, 0),
                    Approx_point(10, 3),
                },
                bbox),
  };

  for(size_t i = 0; i < cases.size(); ++i) {
    const auto& test_case = cases[i];
    test_geom_fill_corners_output_iterator(test_case);
  }
}