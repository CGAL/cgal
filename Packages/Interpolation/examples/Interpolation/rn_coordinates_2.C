//
//file: examples/Interpolation/rn_coordinates_2.C 
//
#include <CGAL/basic.h>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_neighbor_coordinates_traits_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Regular_neighbor_coordinates_traits_2<K>     Gt;
typedef CGAL::Regular_triangulation_2<Gt>     Regular_triangulation;
typedef Regular_triangulation::Weighted_point Weighted_point;

int main()
{
  Regular_triangulation rt;
  
  for (int y=0 ; y<3 ; y++)
    for (int x=0 ; x<3 ; x++) 
      rt.insert(Weighted_point(K::Point_2(x,y), 0));
   
  //coordinate computation
  Weighted_point wp(K::Point_2(1.2, 0.7),2);
  
  std::vector< std::pair< K::Point_2, K::FT  > > coords;
  K::FT  norm = 
    CGAL::regular_neighbor_coordinates_2(rt, wp,
					 std::back_inserter(coords)).second;
  
  return 0; 
}
