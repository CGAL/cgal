#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_polyline_segment_primitive_2.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

typedef std::vector<Point> PointRange;
typedef PointRange::const_iterator Iterator_pr;
typedef CGAL::AABB_polyline_segment_primitive_2<K, Iterator_pr, PointRange> Primitive_pr;
typedef CGAL::AABB_traits_2<K, Primitive_pr> Traits_pr;
typedef CGAL::AABB_tree<Traits_pr> Tree_pr;
typedef Tree_pr::Point_and_primitive_id Point_and_primitive_id_pr;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef Polygon_2::const_iterator Iterator_poly;
typedef CGAL::AABB_polyline_segment_primitive_2<K, Iterator_poly, Polygon_2> Primitive_poly;
typedef CGAL::AABB_traits_2<K, Primitive_poly> Traits_poly;
typedef CGAL::AABB_tree<Traits_poly> Tree_poly;
typedef Tree_poly::Point_and_primitive_id Point_and_primitive_id_poly;

template<class AABBTree, class PPId>
void test(AABBTree tree) {
  tree.build();

  tree.accelerate_distance_queries();

  // counts #intersections with a segment query
  Segment segment_query(Point(1.0, 0.0), Point(0.0, 7.0));

  assert(tree.number_of_intersected_primitives(segment_query) == 2);

  // computes the closest point from a point query
  Point point_query(4.0, 5.0);
  Point closest = tree.closest_point(point_query);
  assert(closest == Point(3.0, 4.0));

  PPId id = tree.closest_point_and_primitive(point_query);
  assert(id.first == closest);
}

int main()
{
  Point a(0.0, 0.0);
  Point b(2.0, 1.0);
  Point c(3.0, 4.0);
  Point d(1.0, 6.0);
  Point e(-1.0, 3.0);

  std::vector<Point> polyline = { a, b, c, d, e };

  Polygon_2 poly(polyline.begin(), polyline.end());

  test<Tree_poly, Point_and_primitive_id_poly>(Tree_poly(poly.begin(), poly.end(), poly));

  test<Tree_pr, Point_and_primitive_id_pr>(Tree_pr(polyline.begin(), std::prev(polyline.end()), polyline));

  return EXIT_SUCCESS;
}
