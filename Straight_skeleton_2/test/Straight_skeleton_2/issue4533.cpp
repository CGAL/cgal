#include<boost/shared_ptr.hpp>
#include <cassert>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_straight_skeleton_2.h>

#include "print.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_2 poly ;


    poly.push_back( Point( 0, 0 ) );
    poly.push_back( Point( 2000, 8000 ) );
    poly.push_back( Point( 10000, 10000 ) );
    poly.push_back( Point( 2000, 12000 ) );
    poly.push_back( Point( 0, 20000 ) );
    poly.push_back( Point( -2000, 12000 ) );
    poly.push_back( Point( -10000, 10000 ) );
    poly.push_back( Point( -2000, 8000 ) );
    assert(poly.is_simple());
    assert(poly.is_counterclockwise_oriented());

  // You can pass the polygon via an iterator pair
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

  print_straight_skeleton(*iss);

  return 0;
}
