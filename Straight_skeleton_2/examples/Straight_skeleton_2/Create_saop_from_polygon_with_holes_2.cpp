#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include "print.h"

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Ss ;

typedef boost::shared_ptr<Polygon_2> PolygonPtr ;
typedef boost::shared_ptr<Ss> SsPtr ;

typedef std::vector<PolygonPtr> PolygonPtrVector ;

int main()
{
  Polygon_2 outer ;
  outer.push_back( Point(-1,-1) ) ;
  outer.push_back( Point(0,-12) ) ;
  outer.push_back( Point(1,-1) ) ;
  outer.push_back( Point(12,0) ) ;
  outer.push_back( Point(1,1) ) ;
  outer.push_back( Point(0,12) ) ;
  outer.push_back( Point(-1,1) ) ;
  outer.push_back( Point(-12,0) ) ;

  Polygon_2 hole ;
  hole.push_back( Point(-1,0) ) ;
  hole.push_back( Point(0,1 ) ) ;
  hole.push_back( Point(1,0 ) ) ;
  hole.push_back( Point(0,-1) ) ;

  assert(outer.is_counterclockwise_oriented());
  assert(hole.is_clockwise_oriented());

  Polygon_with_holes poly( outer ) ;
  poly.add_hole( hole ) ;

  double lOffset = 0.2 ;
  PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset,poly);
  print_polygons(offset_polygons);

  return 0;
}
