
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  std::cout << octree << std::endl;

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  octree.grade();

  std::cout << octree << std::endl;

  std::cout << "\n\n\n" << std::endl;
}

void test_8_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  std::cout << octree << std::endl;

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  octree.grade();

  std::cout << octree << std::endl;

  std::cout << "\n\n\n" << std::endl;
}

void test_10_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});
  points.insert({0.875, 1, -1});
  points.insert({-1, -0.75, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  std::cout << octree << std::endl;

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  octree.grade();

  std::cout << octree << std::endl;

  std::cout << "\n\n\n" << std::endl;
}

int main(void) {

  test_1_point();
  test_8_points();
  test_10_points();

  return EXIT_SUCCESS;
}
