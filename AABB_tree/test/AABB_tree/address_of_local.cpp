#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::Segment_2 Segment_2;
typedef K::Point_2 Point_2;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef Polygon_2::Edge_const_iterator Edge_const_iterator;
typedef CGAL::internal::Source_of_segment_2_iterator_property_map<K,Edge_const_iterator> Sosi_pmap;

int main()
{
  Polygon_2 poly;
  poly.push_back(Point_2(0,0));
  poly.push_back(Point_2(1,0));
  poly.push_back(Point_2(1,1));

  Edge_const_iterator it =      poly.edges_begin();

  Sosi_pmap sosi_pmap;

  Point_2 p = get(sosi_pmap, it);

  return 0;
}
