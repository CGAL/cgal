#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include "dump_to_eps.h"
#include "read_polygon.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Ss ;


typedef boost::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_with_holes input ;

  std::string name = "large_1.poly"; 
  
  std::ifstream is(name.c_str()) ;
  if ( is )
  {
    read_polygon_with_holes(is, input) ;

    boost::optional<double> lMaxTime(2.5);
  
    SsPtr ss = CGAL::create_interior_straight_skeleton_2(input,lMaxTime);

    dump_ss_to_eps(input,ss,"partial_skeleton.eps");
  }
  else
  {
      std::cerr << "Could not open input file: " << name << std::endl;
  }
   
  return 0;
}
