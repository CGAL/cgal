#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include "print.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;


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
    
  Polygon_with_holes poly( outer ) ;
  
  poly.add_hole( hole ) ;
     
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);
  
  print_straight_skeleton(*iss);
  
  return 0;
}
