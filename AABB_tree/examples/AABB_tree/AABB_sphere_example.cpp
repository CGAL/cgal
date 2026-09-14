// Author : Pierre Alliez

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_sphere_primitive_3.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Triangle_3 Triangle;
typedef K::Sphere_3 Sphere;

typedef std::vector<Sphere>::iterator Iterator;
typedef CGAL::AABB_sphere_primitive_3<K, Iterator> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main()
{
  Sphere a(Point(0, 0, 0), 1);
  Sphere b(Point(2, 0, 0), 1.1);
  Sphere c(Point(2, 0, 0), 0.5);

  std::vector<Sphere> spheres = {a, b, c};

  // constructs the AABB tree and the internal search tree for
  // efficient distance computations.
  Tree tree(spheres.begin(), spheres.end());

  // counts #intersections with a plane query
  Plane plane_query(Point(0, 1.04, 0), Point(0, 1.04, 1), Point(1, 1.04, 0));
  std::cout << tree.number_of_intersected_primitives(plane_query)
    << " intersections(s) with plane" << std::endl;

  // counts #intersections with a triangle query
  Triangle triangle_query(Point(0.95, 0, 0), Point(0.95, 1, 1), Point(0.95, 1, 0));
  std::cout << tree.number_of_intersected_primitives(triangle_query)
    << " intersections(s) with triangle" << std::endl;

  // computes the closest point from a point query
  Point point_query(2.0, 2.0, 2.0);
  Tree::Point_and_primitive_id closest = tree.closest_point_and_primitive(point_query);

  std::cout << typeid(*closest.second).name() << std::endl;

  std::cerr << "closest point is: " << closest.first << std::endl;
  return EXIT_SUCCESS;
}
