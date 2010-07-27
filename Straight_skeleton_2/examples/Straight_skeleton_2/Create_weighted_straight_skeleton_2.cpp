#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_straight_skeleton_2.h>

#include "dump_to_eps.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon poly ;
  
  poly.push_back( Point(-1,-1) ) ;
  poly.push_back( Point(0,-12) ) ;
  poly.push_back( Point(1,-1) ) ;
  poly.push_back( Point(12,0) ) ;
  poly.push_back( Point(1,1) ) ;
  poly.push_back( Point(0,12) ) ;
  poly.push_back( Point(-1,1) ) ;
  poly.push_back( Point(-12,0) ) ;
  
  SsPtr ss0 = CGAL::create_straight_skeleton_2(poly, -1.5);

  dump_ss_to_eps(poly,ss0,"uniform_weights_skeleton.eps");

  double weights[] = { -1.0, -1.0, -1.5, -1.0, -0.5, 1.0, 1.5, -1.0 } ;
  
  SsPtr ss1 = CGAL::create_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end(), weights, weights+8);

  dump_ss_to_eps(poly,ss1,"mixed_weights_skeleton.eps");

  return 0;
}
