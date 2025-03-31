#define CGAL_SLS_DEBUG_DRAW

#include <iostream>
#include <iomanip>
#include <string>

#define CGAL_SLS_PRINT_QUEUE_BEFORE_EACH_POP
#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 10000000
#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE 10000000
#define CGAL_STRAIGHT_SKELETON_VALIDITY_ENABLE_TRACE
#define CGAL_POLYGON_OFFSET_ENABLE_TRACE 10000000

void Straight_skeleton_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}

void Straight_skeleton_traits_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <boost/shared_ptr.hpp>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT                        FT ;
typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;

typedef CGAL::Straight_skeleton_2<K> Ss ;
typedef boost::shared_ptr<Ss>        SsPtr ;

int main()
{
#if 0
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
  CGAL::draw(poly);

  // You can pass the polygon via an iterator pair
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
  CGAL::draw(*iss);

  // Or you can pass the polygon directly, as below.

  // To create an exterior straight skeleton you need to specify a maximum offset.
  double lMaxOffset = 5 ;
  SsPtr oss = CGAL::create_exterior_straight_skeleton_2(lMaxOffset, poly);

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*oss);
  CGAL::draw(*oss);
#else
  Polygon_2 poly ;
  poly.push_back( Point(0, 0) );
  // poly.push_back( Point(5, 0) );
  poly.push_back( Point(10, 0) );
  poly.push_back( Point(10, 10) );
  poly.push_back( Point(0, 10) );

  // You can also use weights to modify the speed of some of the fronts
  std::vector<FT> weights = {{ 3, 5, 6, 7 }};
  SsPtr iss = CGAL::create_interior_weighted_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end(),
                                                                 weights.begin(), weights.end());
  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
  CGAL::draw(*iss);
#endif

  return EXIT_SUCCESS;
}
