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

typedef std::list<Segment> SegmentRange;
typedef SegmentRange::iterator Iterator_seg;
typedef CGAL::AABB_segment_primitive_2<K, Iterator_seg> Primitive_seg;
typedef CGAL::AABB_traits_2<K, Primitive_seg> Traits_seg;
typedef CGAL::AABB_tree<Traits_seg> Tree_seg;

typedef std::vector<Point> PointRange;
typedef PointRange::iterator Iterator_pr;
typedef CGAL::AABB_polyline_segment_primitive_2<PointRange, K, Iterator_pr> Primitive_pr;
typedef CGAL::AABB_traits_2<K, Primitive_pr> Traits_pr;
typedef CGAL::AABB_tree<Traits_pr> Tree_pr;

typedef CGAL::Polygon_2<K> Polygon;
typedef Polygon::iterator Iterator_poly;
typedef CGAL::AABB_polyline_segment_primitive_2<Polygon, K, Iterator_poly> Primitive_poly;
typedef CGAL::AABB_traits_2<K, Primitive_poly> Traits_poly;
typedef CGAL::AABB_tree<Traits_poly> Tree_poly;


int main()
{
  Point a(0.0, 0.0);
  Point b(2.0, 1.0);
  Point c(3.0, 4.0);
  Point d(1.0, 6.0);
  Point e(-1.0, 3.0);

  std::list<Segment> seg;
  seg.push_back(Segment(a, b));
  seg.push_back(Segment(b, c));
  seg.push_back(Segment(c, d));
  seg.push_back(Segment(d, e));
  seg.push_back(Segment(e, a));

  std::vector<Point> polyline = { a, b, c, d, e };

  Polygon poly(polyline.begin(), polyline.end());

  //   // constructs the AABB tree and the internal search tree for
  //   // efficient distance computations.
  Tree_seg tree_seg(seg.begin(), seg.end());
  Tree_pr tree_pr(polyline.begin(), polyline.end(), polyline);
  Tree_poly tree_poly(poly.begin(), poly.end(), poly);

  tree_seg.accelerate_distance_queries();
  tree_pr.accelerate_distance_queries();
  tree_poly.accelerate_distance_queries();

  //   // counts #intersections with a segment query

  Segment segment_query(Point(1.0, 0.0), Point(0.0, 7.0));
  std::cout << tree_seg.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  std::cout << tree_pr.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  std::cout << tree_poly.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  // computes the closest point from a point query
  Point point_query(1.5, 3.0);
  Point closest = tree_seg.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;
  closest = tree_pr.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;
  closest = tree_poly.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;
  return EXIT_SUCCESS;
}
