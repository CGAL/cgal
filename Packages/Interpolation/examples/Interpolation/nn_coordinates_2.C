//
//file: examples/Interpolation/nn_coordinates_2.C 
//
#include <CGAL/basic.h>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;

int main()
{
  Delaunay_triangulation dt;
  
  for (int y=0 ; y<3 ; y++)
    for (int x=0 ; x<3 ; x++) 
      dt.insert(K::Point_2(x,y));
   
  //coordinate computation
  K::Point_2 p(1.2, 0.7);
  std::vector< std::pair< K::Point_2, K::FT  > > coords;
  K::FT  norm = 
    CGAL::natural_neighbor_coordinates_2(dt, p,
					 std::back_inserter(coords)).second;
  
  return 0; 
}
