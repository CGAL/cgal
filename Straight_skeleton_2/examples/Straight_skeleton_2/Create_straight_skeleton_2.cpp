#include<vector>
#include<iterator>
#include<iostream>
#include<iomanip>
#include<string>

#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/Create_straight_skeleton_2.h>

#include "print.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point_2 ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_2 poly ;
  
  poly.push_back( Point_2(-1,-1) ) ;
  poly.push_back( Point_2(0,-12) ) ;
  poly.push_back( Point_2(1,-1) ) ;
  poly.push_back( Point_2(12,0) ) ;
  poly.push_back( Point_2(1,1) ) ;
  poly.push_back( Point_2(0,12) ) ;
  poly.push_back( Point_2(-1,1) ) ;
  poly.push_back( Point_2(-12,0) ) ;
     
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);

  print_straight_skeleton(*iss);
  
  // To create an exterior straight skeleton you need to specify a maximum offset
  double lMaxOffset = 5 ; 
  SsPtr oss = CGAL::create_exterior_straight_skeleton_2(lMaxOffset,poly);
  
  print_straight_skeleton(*oss);
  
  return 0;
}
