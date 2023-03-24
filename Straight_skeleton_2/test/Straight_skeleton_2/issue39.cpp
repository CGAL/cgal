#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print_straight_skeleton_2.h>
#include <CGAL/draw_straight_skeleton_2.h>

#include <boost/shared_ptr.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

int main(int, char**)
{
  Polygon_2 poly ;

  poly.push_back( Point(320, 1120) ) ;
  poly.push_back( Point(-141.55499021951055738, 1658.4808219227620611) ) ;
  poly.push_back( Point(-1115.9820207544369168, 1835.6551586833445526) ) ;
  poly.push_back( Point(760, 520) ) ;
  poly.push_back( Point(800, 560) ) ;
  poly.push_back( Point(840, 600) ) ;
  poly.push_back( Point(800,640) ) ;
  poly.push_back( Point(400, 1040) ) ;

  // You can pass the polygon via an iterator pair
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());
  if( !iss )
  {
    std::cout << "Failed to create interior straight skeleton" << std::endl ;
    return EXIT_FAILURE ;
  }

  // Or you can pass the polygon directly, as below.

  // To create an exterior straight skeleton you need to specify a maximum offset.
  double lMaxOffset = 5 ;
  SsPtr oss = CGAL::create_exterior_straight_skeleton_2(lMaxOffset, poly);

  if( !oss )
  {
    std::cout << "Failed to create exterior straight skeleton" << std::endl ;
    return EXIT_FAILURE ;
  }

  print_straight_skeleton(*iss);
  print_straight_skeleton(*oss);

  CGAL::draw(*iss);
  CGAL::draw(*oss);

  std::cout << "OK" << std::endl;

  return EXIT_SUCCESS ;
}
