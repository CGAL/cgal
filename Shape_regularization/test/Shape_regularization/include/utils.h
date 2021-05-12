#ifndef CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H

// STL includes.
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/Random.h>

namespace CGAL {
namespace Shape_regularization {
namespace Tests {

template<typename Segment_2>
void check_reference_values(
  const std::vector<Segment_2>& input_range,
  const std::vector<int>& reference_values) {

  for (std::size_t i = 0; i < input_range.size(); ++i) {
    const auto& segment = input_range[i];

    const auto point1 = segment.source().x() + segment.source().y();
    const auto point2 = segment.target().x() + segment.target().y();

    const int key = static_cast<int>(floor(
      CGAL::to_double(point1 + point2)));
    // std::cout << key << std::endl;
    assert(key == reference_values[i]);
  }
}

double get_coefficient_value(
  const double theta, double& iterator) {

  if (
    theta == 0.0 ||
    theta == CGAL_PI / 2.0 ||
    theta == CGAL_PI ||
    theta == 3.0 * CGAL_PI / 2.0) {

    iterator = 0.0;
  } else if (
    theta == CGAL_PI / 4.0 ||
    theta == 3.0 * CGAL_PI / 4.0 ||
    theta == 5.0 * CGAL_PI / 4.0 ||
    theta == 7.0 * CGAL_PI / 4.0) {

    iterator = 0.22;
  } else if (
    (theta > 0.0 && theta < CGAL_PI / 4.0) ||
    (theta > CGAL_PI / 2.0 && theta < 3.0 * CGAL_PI / 4.0) ||
    (theta > CGAL_PI && theta < 5.0 * CGAL_PI / 4.0) ||
    (theta > 3.0 * CGAL_PI / 2.0 && theta < 7.0 * CGAL_PI / 4.0)) {

    iterator += 0.02;
  } else
    iterator -= 0.02;

  if (theta < CGAL_PI) return -1.0 * iterator;
  return iterator;
}

template<typename Segment_2>
void create_example_offsets(
  std::vector<Segment_2>& segments) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using Point_2 = typename Traits::Point_2;

  segments.clear();
  segments.reserve(100);

  double theta = 0.0, coef = 0.0, iterator = 0.0;
  double theta_step = CGAL_PI / 25.0;

  while (theta < 2.0 * CGAL_PI) {
    const double st = std::sin(theta);
    const double ct = std::cos(theta);

    const Point_2 a = Point_2(0.0, 0.0);
    const Point_2 b = Point_2(ct, st);

    coef = get_coefficient_value(theta, iterator);
    const Point_2 c = Point_2(ct, st + coef);
    const Point_2 d = Point_2(2.0 * ct, 2.0 * st + coef);
    theta += theta_step;

    segments.push_back(Segment_2(a, b));
    segments.push_back(Segment_2(c, d));
  }
}

template<typename Segment_2>
void create_example_angles(
  std::vector<Segment_2>& segments) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using Point_2 = typename Traits::Point_2;

  using PG = CGAL::Points_on_segment_2<Point_2>;
  using Creator = CGAL::Creator_uniform_2<Point_2, Segment_2>;
  using Segment_iterator = CGAL::Join_input_iterator_2<PG, PG, Creator>;
  using Count_iterator = CGAL::Counting_iterator<Segment_iterator, Segment_2>;

  segments.clear();
  segments.reserve(100);

  // A horizontal like fan.
  PG p1(Point_2(-250,  -50), Point_2(-250,  50), 50);
  PG p2(Point_2( 250, -250), Point_2( 250, 250), 50);

  Segment_iterator t1(p1, p2);
  Count_iterator t1_begin(t1);
  Count_iterator t1_end(t1, 50);
  std::copy(t1_begin, t1_end, std::back_inserter(segments));

  // A vertical like fan.
  PG p3(Point_2( -50, -250), Point_2( 50, -250), 50);
  PG p4(Point_2(-250,  250), Point_2(250,  250), 50);

  Segment_iterator t2(p3, p4);
  Count_iterator t2_begin(t2);
  Count_iterator t2_end(t2, 50);
  std::copy(t2_begin, t2_end, std::back_inserter(segments));
}

template<typename Point_2>
void create_example_multiple_directions(
  std::vector<Point_2>& contour) {

  contour.clear();
  contour = {
    Point_2(546115.55975486419629, 5244326.4265436409041),
    Point_2(546117.37262898986228, 5244337.1341805141419),
    Point_2(546116.01587088918313, 5244341.0215584803373),
    Point_2(546117.94300083850976, 5244350.5947591038421),
    Point_2(546112.84907653776463, 5244365.2875102143735),
    Point_2(546107.19443975773174, 5244363.9731347402558),
    Point_2(546106.22850909293629, 5244364.9865339472890),
    Point_2(546110.05300284817349, 5244372.5691289855167),
    Point_2(546105.07976105506532, 5244374.5981351630762),
    Point_2(546093.72503010707442, 5244352.3779847342521),
    Point_2(546091.96853454655502, 5244351.3945095967501),
    Point_2(546091.50930348574184, 5244352.8908741353080),
    Point_2(546103.43633847369347, 5244375.3802235340700),
    Point_2(546102.20120883965865, 5244376.9659267812967),
    Point_2(546091.04900722880848, 5244382.3613784974441),
    Point_2(546081.03281952394173, 5244362.8676363257691),
    Point_2(546078.12493466213346, 5244361.0630018385127),
    Point_2(546078.89294813224114, 5244359.1702114883810),
    Point_2(546073.78569827671163, 5244360.0353698069230),
    Point_2(546066.88938752864487, 5244346.7618661979213),
    Point_2(546088.89068568043876, 5244335.1294597545639),
    Point_2(546094.60000572062563, 5244346.0973066119477),
    Point_2(546096.67206490959506, 5244346.9671173477545),
    Point_2(546103.79629588325042, 5244359.4217483354732),
    Point_2(546107.84054918668699, 5244348.1492590270936),
    Point_2(546106.88664358214010, 5244340.8161107189953),
    Point_2(546100.48662372061517, 5244340.8635002989322),
    Point_2(546098.61382249533199, 5244329.1317009674385),
    Point_2(546106.99944098747801, 5244326.3489141445607)
  };
}

} // namespace Tests
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H
