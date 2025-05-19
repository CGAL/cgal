#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <memory>

#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT                        FT ;
typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef std::shared_ptr<Polygon_2> PolygonPtr ;
typedef std::shared_ptr<Ss> SsPtr ;

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

  FT lOffset = 1 ;

  PolygonPtrVector inner_offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset,poly);
  PolygonPtrVector outer_offset_polygons = CGAL::create_exterior_skeleton_and_offset_polygons_2(lOffset,poly);

  CGAL::Straight_skeletons_2::IO::print_polygons(inner_offset_polygons);
  CGAL::Straight_skeletons_2::IO::print_polygons(outer_offset_polygons);

  return EXIT_SUCCESS;
}
