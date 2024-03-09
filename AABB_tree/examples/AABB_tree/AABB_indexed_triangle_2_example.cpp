// Author(s) : Camille Wormser, Pierre Alliez

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_indexed_triangle_primitive_2.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_2 Point;

typedef std::vector<std::array<uint8_t, 3> >::iterator IndexIterator;
typedef std::vector<Point> PointRange;
typedef CGAL::AABB_indexed_triangle_primitive_2<K, IndexIterator, PointRange> Primitive;
typedef CGAL::AABB_traits_2<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

int main()
{
  Point a(0.0, 0.0);
  Point b(0.0, 1.0);
  Point c(1.0, 0.0);
  Point d(1.0, 1.0);
  Point e(2.0, 0.0);
  Point f(2.0, 1.0);

  std::vector<Point> points = { a, b, c, d, e, f };

  std::vector<std::array<uint8_t, 3> > triangles;
  triangles.push_back({ 0, 2, 1 });
  triangles.push_back({ 1, 2, 3 });
  triangles.push_back({ 3, 2, 4 });
  triangles.push_back({ 3, 4, 5 });

  // constructs AABB tree
  Tree tree(triangles.begin(), triangles.end(), points);

  // point sampling
  Point_and_primitive_id id;
  id = tree.closest_point_and_primitive(Point(0.5, 0.4));
  std::cout << std::distance(triangles.begin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point(0.5, 0.6));
  std::cout << std::distance(triangles.begin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point(1.5, 0.5));
  std::cout << std::distance(triangles.begin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point(1.5, 0.6));
  std::cout << std::distance(triangles.begin(), id.second) << ". triangle" << std::endl;

  return EXIT_SUCCESS;
}
