#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::Segment_2 Segment_2;
typedef K::Point_2 Point_2;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef Polygon_2::Edge_const_iterator Edge_const_iterator;
typedef CGAL::internal::Source_of_segment_2_iterator_property_map<K,Edge_const_iterator> Sosi_polygon_map;

typedef std::vector<Segment_2>::const_iterator Seg_const_iterator;
typedef CGAL::internal::Source_of_segment_2_iterator_property_map<K,Seg_const_iterator> Sosi_segment_map;

int main()
{
  Polygon_2 poly;
  poly.push_back(Point_2(0,0));
  poly.push_back(Point_2(1,0));
  poly.push_back(Point_2(1,1));

  Edge_const_iterator it =      poly.edges_begin();

  Sosi_polygon_map sosi_polygon_map;

  std::vector<Segment_2> segs;
  segs.push_back(Segment_2(Point_2(0,0), Point_2(1,0)));
  segs.push_back(Segment_2(Point_2(1,0), Point_2(1,1)));
  segs.push_back(Segment_2(Point_2(1,1), Point_2(1,0)));

  Sosi_segment_map sosi_segment_map;

  const Point_2& p = get(sosi_polygon_map, it);
  const Point_2& q = get(sosi_segment_map, segs.begin());
  assert(p == poly[0]);
  assert(q == segs[0].source());
  return 0;
}
