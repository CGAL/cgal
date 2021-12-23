#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include "print.h"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Polygon_2> PolygonPtr ;
typedef boost::shared_ptr<Ss> SsPtr ;

typedef std::vector<PolygonPtr> PolygonPtrVector ;

int main()
{
  Polygon_2 poly ;
  poly.push_back( Point(-1,-1) ) ;
  poly.push_back( Point(0,-12) ) ;
  poly.push_back( Point(1,-1) ) ;
  poly.push_back( Point(12,0) ) ;
  poly.push_back( Point(1,1) ) ;
  poly.push_back( Point(0,12) ) ;
  poly.push_back( Point(-1,1) ) ;
  poly.push_back( Point(-12,0) ) ;
  assert(poly.is_counterclockwise_oriented());

  SsPtr ss = CGAL::create_interior_straight_skeleton_2(poly);

  double lOffset = 1 ;
  PolygonPtrVector offset_polygons = CGAL::create_offset_polygons_2<Polygon_2>(lOffset,*ss);
  print_polygons(offset_polygons);

  return 0;
}
