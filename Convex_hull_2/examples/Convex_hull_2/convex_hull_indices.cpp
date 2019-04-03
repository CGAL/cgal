#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Convex_hull_traits_adapter_2<K,
          CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;

int main()
{
  std::vector<Point_2> points;

  points.push_back(Point_2(10,0));
  points.push_back(Point_2(0,10));
  points.push_back(Point_2(1,1));
  points.push_back(Point_2(3,4));
  points.push_back(Point_2(0,0));
  points.push_back(Point_2(10,10));
  points.push_back(Point_2(2,6));
  
  std::vector<std::size_t> indices(points.size()), out;

  for(int i=0; i < indices.size();i++){
    indices[i] = i;
  }

  
  CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                      Convex_hull_traits_2(CGAL::make_property_map(points)));

  for( std::size_t i : out){
    std::cout << "points[" << i << "] = " << points[i] << std::endl;
  }

  return 0;
}
