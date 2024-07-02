#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

typedef std::list<Segment> SegmentRange;
typedef SegmentRange::const_iterator Iterator;
typedef CGAL::AABB_segment_primitive_2<K, Iterator> Primitive;
typedef CGAL::AABB_traits_2<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

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

  // constructs the AABB tree and the internal search tree for
  // efficient distance computations.
  Tree tree(seg.begin(), seg.end());
  tree.build();

  tree.accelerate_distance_queries();

  // counts #intersections with a segment query
  Segment segment_query(Point(1.0, 0.0), Point(0.0, 7.0));
  std::cout << tree.number_of_intersected_primitives(segment_query)
    << " intersections(s) with segment" << std::endl;

  // computes the closest point from a point query
  Point point_query(1.5, 3.0);
  Point closest = tree.closest_point(point_query);
  std::cerr << "closest point is: " << closest << std::endl;

  Point_and_primitive_id id = tree.closest_point_and_primitive(point_query);
  std::cout << id.second->source() << " " << id.second->target() << std::endl;

  return EXIT_SUCCESS;
}
