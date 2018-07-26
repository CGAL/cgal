#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point_3;

typedef CGAL::Surface_mesh<Point_3>                             Mesh;
typedef std::vector<std::size_t>                                Polygon;

void test_simplify_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test simplify_polygons... " << std::endl;

  std::vector<Point_3> points;
  std::vector<Polygon> polygons;

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,2,0)); // #1
  points.push_back(Point_3(1,0,0)); // #2
  points.push_back(Point_3(1,3,0)); // #3
  points.push_back(Point_3(0,1,0)); // #4
  points.push_back(Point_3(1,1,0)); // #5
  points.push_back(Point_3(0,0,0)); // #6

  // ------
  Polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  std::size_t res = PMP::internal::simplify_polygons_in_polygon_soup<K>(points, polygons);
  assert(res == 0 && polygons.back().size() == 3);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 1);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(0); polygon.push_back(0); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 1);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(3); polygon.push_back(3); polygon.push_back(3); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 2);

  // ------
  // Now with the same geometric positions, but different combinatorial information
  polygon.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(1); polygon.push_back(6);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 3);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(6); polygon.push_back(0); polygon.push_back(5);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 2);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(6); polygon.push_back(5); polygon.push_back(3); polygon.push_back(3);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.back().size() == 3);
}

void test_remove_invalid_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test remove_invalid_polygons... " << std::endl;

  // points are not actually needed since only the size of the polygons is considered
  std::vector<Point_3> points;
  std::vector<Polygon> polygons;

  std::size_t res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 0 && polygons.size() == 0);

  // non-trivial polygon
  Polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 0 && polygons.size() == 1);

  // another non-trivial polygon
  polygon.clear();
  polygon.push_back(1); polygon.push_back(3); polygon.push_back(3); polygon.push_back(7);
  polygons.push_back(polygon);

  // empty polygon
  polygon.clear();
  polygons.push_back(polygon);

  // 1-vertex polygon
  polygon.push_back(0);
  polygons.push_back(polygon);

  // 2-vertex polygon
  polygon.push_back(1);
  polygons.push_back(polygon);

  // another non-trivial polygon
  polygon.clear();
  polygon.push_back(8); polygon.push_back(9); polygon.push_back(7); polygon.push_back(6);
  polygons.push_back(polygon);

  res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 3 && polygons.size() == 3);
}

template <typename PointRange, typename PolygonRange>
std::size_t test_remove_isolated_points_data_set(PointRange& points,
                                                 PolygonRange& polygons,
                                                 const bool verbose = false)
{
  if(verbose)
  {
    std::cout << "Input:" << std::endl;
    std::cout << points.size() << " points" << std::endl;
    std::cout << polygons.size() << " polygons" << std::endl;
  }

  std::size_t rm_nv = PMP::remove_isolated_points_in_polygon_soup(points, polygons);

  if(verbose)
  {
    std::cout << "Removed " << rm_nv << " points" << std::endl;
    std::cout << points.size() << " points" << std::endl;
    std::cout << polygons.size() << " polygons" << std::endl;

    std::cout << "Polygons:" << std::endl;
    for(std::size_t i=0; i<polygons.size(); ++i)
    {
      for(std::size_t j=0; j<polygons[i].size(); ++j)
        std::cout << polygons[i][j] << " ";
      std::cout << std::endl;
    }
  }

  return rm_nv;
}

void test_remove_isolated_points(const bool verbose = false)
{
  std::cout << "test remove_isolated_points... " << std::endl;

  std::vector<Point_3> points;
  std::vector<Polygon> polygons;

  // everything empty
  std::size_t res = test_remove_isolated_points_data_set(points, polygons, verbose);
  assert(res == 0 && points.empty() && polygons.empty());

  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,2,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(1,3,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(1,1,0));

  // no polygons (all points are unused)
  std::deque<Point_3> points_copy(points.begin(), points.end());
  res = test_remove_isolated_points_data_set(points_copy, polygons, verbose);
  assert(res == 6 && points_copy.empty() && polygons.empty());

  Polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  polygon.clear();
  polygon.push_back(4); polygon.push_back(0); polygon.push_back(2); polygon.push_back(5);
  polygons.push_back(polygon);

  // generic test
  res = test_remove_isolated_points_data_set(points, polygons, verbose);
  assert(res == 2 && points.size() == 4 && polygons.size() == 2);
}

void test_slit_pinched_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test split_pinched_polygons... " << std::endl;

  std::vector<Point_3> points;
  std::vector<Polygon> polygons;

  // everything empty
  std::size_t res = PMP::internal::split_pinched_polygons_in_polygon_soup<K>(points, polygons);
  assert(res == 0 && points.empty() && polygons.empty());

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,0,0)); // #1
  points.push_back(Point_3(0,1,0)); // #2
  points.push_back(Point_3(1,1,0)); // #3
  points.push_back(Point_3(1,3,0)); // #4
  points.push_back(Point_3(1,1,0)); // #5

  // no pinch
  Polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 0 && polygons.size() == 1);

  // pinch via same ID (2)
  polygon.clear();
  polygons.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4); polygon.push_back(3); polygon.push_back(2);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 3 && polygons[1].size() == 2);

  // pinch via same point (5 & 3)
  polygon.clear();
  polygons.clear();
  polygon.push_back(5); polygon.push_back(1); polygon.push_back(0);
  polygon.push_back(3); polygon.push_back(4); polygon.push_back(2);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 3 && polygons[1].size() == 3);

  // pinch on last point
  polygon.clear();
  polygons.clear();
  polygon.push_back(1); polygon.push_back(5); polygon.push_back(0);
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(3);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 4 && polygons[1].size() == 2);

  // multiple pinches
  polygon.clear();
  polygons.clear();
  polygon.push_back(5); polygon.push_back(1); polygon.push_back(3); // pinch 5&3, pinch 3...
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(0);
  polygon.push_back(1); polygon.push_back(3); polygon.push_back(1); // ... and 3, pinch 1...
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(1); // ... and 1
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 3 && polygons.size() == 4);
  assert(polygons[0].size() == 2); // 5 1
  assert(polygons[1].size() == 5); // 3 2 4 5 1
  assert(polygons[2].size() == 3); // 1 2 4
  assert(polygons[3].size() == 2); // 1 3
}

int main()
{
  test_simplify_polygons(false);
  test_remove_invalid_polygons(false);
  test_remove_isolated_points(false);
  test_slit_pinched_polygons(false);

  return EXIT_SUCCESS;
}
