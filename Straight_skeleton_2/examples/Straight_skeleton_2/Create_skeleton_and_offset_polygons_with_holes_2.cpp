#include<vector>

#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_offset_polygons_from_polygon_with_holes_2.h>

#include "print.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> PolygonWithHoles ;

typedef boost::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr ;

typedef std::vector<PolygonWithHolesPtr> PolygonWithHolesPtrVector;


int main()
{  
  Polygon_2 outer ;  
  
  outer.push_back( Point( 0.0, 0.0) ) ;
  outer.push_back( Point(10.0, 0.0) ) ;
  outer.push_back( Point(10.0, 4.5) ) ;
  outer.push_back( Point(12.0, 4.5) ) ;
  outer.push_back( Point(12.0, 2.0) ) ;
  outer.push_back( Point(16.0, 2.0) ) ;
  outer.push_back( Point(16.0, 8.0) ) ;
  outer.push_back( Point(12.0, 8.0) ) ;
  outer.push_back( Point(12.0, 5.5) ) ;
  outer.push_back( Point(10.0, 5.5) ) ;
  outer.push_back( Point(10.0,10.0) ) ;
  outer.push_back( Point( 0.0,10.0) ) ;
  
  Polygon_2 hole ;
  
  hole.push_back( Point(3.0,3.0) ) ;
  hole.push_back( Point(3.0,7.0) ) ;
  hole.push_back( Point(7.0,7.0) ) ;
  hole.push_back( Point(7.0,3.0) ) ;
    
  PolygonWithHoles poly( outer ) ;
  
  poly.add_hole( hole ) ;
      
  double lOffset = 1 ;

  PolygonWithHolesPtrVector offset_poly_with_holes = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(lOffset,poly);
  
  print_polygons_with_holes(offset_poly_with_holes);

  return 0;
}
