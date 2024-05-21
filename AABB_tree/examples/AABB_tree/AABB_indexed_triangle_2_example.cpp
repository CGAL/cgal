#include <iostream>
#include <vector>
#include <array>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_indexed_triangle_primitive_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef K::Ray_2 Ray;

template <class GeomTraits>
struct Projection_xy_point_map {

  typedef typename GeomTraits::Point_3 key_type;
  typedef typename GeomTraits::Point_2 value_type;
  typedef value_type reference;

  typedef boost::readable_property_map_tag category;
  typedef Projection_xy_point_map<GeomTraits> Self;

  Projection_xy_point_map() {}

  inline friend value_type get(Self, key_type p)
  {
    return value_type(p.x(), p.y());
  }
};

typedef std::vector<std::array<uint8_t, 3> >::const_iterator IndexIterator;
typedef std::vector<Point_3> PointRange;
typedef CGAL::AABB_indexed_triangle_primitive_2<K, IndexIterator, PointRange, CGAL::Tag_false, Projection_xy_point_map<K>> Primitive;
typedef CGAL::AABB_traits_2<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef std::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

int main()
{
  Point_3 a(0.0, 0.0, 0.0);
  Point_3 b(0.0, 1.0, 0.0);
  Point_3 c(1.0, 0.0, 0.0);
  Point_3 d(1.0, 1.0, 0.0);
  Point_3 e(2.0, 0.0, 0.0);
  Point_3 f(2.0, 1.0, 0.0);

  std::vector<Point_3> points = { a, b, c, d, e, f };

  std::vector<std::array<uint8_t, 3> > triangles;
  triangles.push_back({ 0, 2, 1 });
  triangles.push_back({ 1, 2, 3 });
  triangles.push_back({ 3, 2, 4 });
  triangles.push_back({ 3, 4, 5 });

  // constructs AABB tree
  Tree tree(triangles.begin(), triangles.end(), points);

  // point sampling
  Point_and_primitive_id id;
  id = tree.closest_point_and_primitive(Point_2(0.5, 0.4));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point_2(0.5, 0.6));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point_2(1.5, 0.5));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point_2(1.5, 0.6));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point_2(1.0, 0.0));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;
  id = tree.closest_point_and_primitive(Point_2(3.0, 0.5));
  std::cout << std::distance(triangles.cbegin(), id.second) << ". triangle" << std::endl;

  Ray ray(Point_2(5.5, 0.5), Point_2(1.5, 0.5));
  Ray_intersection intersection = tree.first_intersection(ray);

  if (!intersection) {
    std::cout << "Ray does not intersect with triangles although it should!" << std::endl;
    return EXIT_FAILURE;
  }
  else {
    std::cout << std::distance(triangles.cbegin(), intersection->second) << ". triangle" << std::endl;
  }

  std::list<Ray_intersection> intersections;
  tree.all_intersections(ray, std::back_inserter(intersections));

  return EXIT_SUCCESS;
}
