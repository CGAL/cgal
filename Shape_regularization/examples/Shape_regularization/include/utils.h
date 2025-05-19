#ifndef CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H

// STL includes.
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/Random.h>
#include <CGAL/IO/io.h>

template<typename Point_2>
void initialize_contour(
  const std::string path,
  std::vector<Point_2>& contour) {

  using Traits = typename CGAL::Kernel_traits<Point_2>::Kernel;
  using FT = typename Traits::FT;

  contour.clear();
  std::ifstream file(path.c_str(), std::ios_base::in);
  CGAL::IO::set_ascii_mode(file);
  file.precision(20);

  if (!file) {
    std::cerr <<
    "Error: cannot read the file with data!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  double xd, yd;
  std::string tmp;
  while (!file.eof()) {
    file >> tmp >> xd >> yd >> tmp >> tmp >> tmp >> tmp;
    const FT x = static_cast<FT>(xd);
    const FT y = static_cast<FT>(yd);
    const Point_2 point = Point_2(x, y);
    contour.push_back(point);
  }
  contour.pop_back();
  file.close();
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
  } else {
    iterator -= 0.02;
  }

  if (theta < CGAL_PI) return -1.0 * iterator;
  return iterator;
}

template<typename Segment_2>
void create_example_offsets(
  std::vector<Segment_2>& segments) {

  using Traits = typename CGAL::Kernel_traits<Segment_2>::Kernel;
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

  using Traits = typename CGAL::Kernel_traits<Segment_2>::Kernel;
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

template<typename Segment_2>
void create_example_15(
  std::vector<Segment_2>& segments) {

  using Traits = typename CGAL::Kernel_traits<Segment_2>::Kernel;
  using Point_2 = typename Traits::Point_2;

  const std::vector<Point_2> points = {
    Point_2(1.000000, 1.000000), Point_2(0.925377, 2.995179),
    Point_2(1.000000, 3.000000), Point_2(1.066662, 4.951894),
    Point_2(1.000000, 5.000000), Point_2(2.950000, 4.930389),
    Point_2(3.000000, 4.950000), Point_2(2.934996, 3.008203),
    Point_2(3.085452, 3.003266), Point_2(2.969782, 1.002004),
    Point_2(0.948866, 3.033161), Point_2(2.900000, 3.000000),
    Point_2(0.930000, 1.000000), Point_2(2.860000, 1.002004),
    Point_2(1.600000, 4.000000), Point_2(1.932136, 4.364718),
    Point_2(1.598613, 3.982686), Point_2(2.018220, 3.686595),
    Point_2(1.951872, 4.363094), Point_2(2.290848, 4.054154),
    Point_2(2.018220, 3.686595), Point_2(2.304517, 4.045054),
    Point_2(1.642059, 1.928505), Point_2(1.993860, 2.247986),
    Point_2(1.993860, 2.247986), Point_2(2.259099, 1.919966),
    Point_2(1.629845, 1.923077), Point_2(1.968759, 1.599174),
    Point_2(2.259099, 1.919966), Point_2(1.968759, 1.599174)
  };

  segments.clear();
  segments.reserve(15);

  segments.push_back(Segment_2(points[0] ,  points[1]));
  segments.push_back(Segment_2(points[2] ,  points[3]));
  segments.push_back(Segment_2(points[4] ,  points[5]));
  segments.push_back(Segment_2(points[6] ,  points[7]));
  segments.push_back(Segment_2(points[8] ,  points[9]));
  segments.push_back(Segment_2(points[10], points[11]));
  segments.push_back(Segment_2(points[12], points[13]));
  segments.push_back(Segment_2(points[14], points[15]));
  segments.push_back(Segment_2(points[16], points[17]));
  segments.push_back(Segment_2(points[18], points[19]));
  segments.push_back(Segment_2(points[20], points[21]));
  segments.push_back(Segment_2(points[22], points[23]));
  segments.push_back(Segment_2(points[24], points[25]));
  segments.push_back(Segment_2(points[26], points[27]));
  segments.push_back(Segment_2(points[28], points[29]));
}

template<
typename Point_2,
typename Line_2>
void boundary_points_on_line_2(
  const std::vector<Point_2>& points,
  const Line_2& line,
  Point_2& p, Point_2& q) {

  using Traits = typename CGAL::Kernel_traits<Point_2>::Kernel;
  using FT = typename Traits::FT;
  using Vector_2 = typename Traits::Vector_2;

  const FT max_value = FT(1000000000000);

  FT min_proj_value =  max_value;
  FT max_proj_value = -max_value;

  const Vector_2 ref_vector = line.to_vector();
  const Point_2& ref_point  = points[0];

  for (const auto& point : points) {
    const Vector_2 curr_vector(ref_point, point);
    const FT value = CGAL::scalar_product(curr_vector, ref_vector);
    if (value < min_proj_value) {
      min_proj_value = value;
      p = point; }
    if (value > max_proj_value) {
      max_proj_value = value;
      q = point; }
  }
}

template<typename Point_2>
void initialize_groups(
  const std::string path,
  std::vector< std::vector<Point_2> >& groups) {

  std::ifstream file(path.c_str(), std::ios_base::in);
  CGAL::IO::set_ascii_mode(file);
  file.precision(20);

  if (!file) {
    std::cout <<
    "Error: cannot read the file with data!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  Point_2 p; double stub;
  std::size_t group_index;
  std::map<std::size_t, std::vector<Point_2> > data_map;
  while (!file.eof()) {
    file >> p >> stub >> group_index;
    data_map[group_index].push_back(p);
  }
  file.close();

  groups.clear();
  groups.reserve(data_map.size());
  for (const auto& pair : data_map)
    groups.push_back(pair.second);
}

#endif // CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H
