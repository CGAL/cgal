#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_polyline_segment_primitive_2.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/Polygon_2.h>


typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

typedef std::vector<Point> PointRange;
typedef PointRange::iterator Iterator_pr;
typedef CGAL::AABB_polyline_segment_primitive_2<K, Iterator_pr, PointRange> Primitive_pr;
typedef CGAL::AABB_traits_2<K, Primitive_pr> Traits_pr;
typedef CGAL::AABB_tree<Traits_pr> Tree_pr;
typedef Tree_pr::Point_and_primitive_id Point_and_primitive_id_pr;

typedef CGAL::Polygon_2<K> Polygon;
typedef Polygon::iterator Iterator_poly;
typedef CGAL::AABB_polyline_segment_primitive_2<K, Iterator_poly, Polygon> Primitive_poly;
typedef CGAL::AABB_traits_2<K, Primitive_poly> Traits_poly;
typedef CGAL::AABB_tree<Traits_poly> Tree_poly;
typedef Tree_poly::Point_and_primitive_id Point_and_primitive_id_poly;


int main()
{
  Point a(0.0, 0.0);
  Point b(2.0, 1.0);
  Point c(3.0, 4.0);
  Point d(1.0, 6.0);
  Point e(-1.0, 3.0);

  std::vector<Point> polyline = { a, b, c, d, e };

  Polygon poly(polyline.begin(), polyline.end());

  // constructs the AABB tree and the internal search tree for
  // efficient distance computations.

  // For a point range, the second iterator must be of the second-last point in the range.
  // If it points to the end of the range, to polyline is considered to be closed.
  Tree_pr tree_pr(polyline.begin(), std::prev(polyline.end()), polyline);
  Tree_poly tree_poly(poly.begin(), poly.end(), poly);

  tree_pr.build();
  tree_poly.build();

  tree_pr.accelerate_distance_queries();
  tree_poly.accelerate_distance_queries();

  // counts #intersections with a segment query
  Segment segment_query(Point(1.0, 0.0), Point(0.0, 7.0));

  std::cout << tree_pr.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  std::cout << tree_poly.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  // computes the closest point from a point query
  Point point_query(1.5, 3.0);
  Point closest = tree_pr.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;

  closest = tree_poly.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;

  Point_and_primitive_id_poly id_poly = tree_poly.closest_point_and_primitive(point_query);
  // id_poly.second is of type Iterator_poly and points to the first point of the segment in the polygon point range.

  Point_and_primitive_id_pr id_pr = tree_pr.closest_point_and_primitive(point_query);
  // id_pr.second is of type Iterator_poly and points to the first point of the segment in the point range.

  return EXIT_SUCCESS;
}
