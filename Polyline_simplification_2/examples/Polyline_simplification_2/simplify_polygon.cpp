#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                   Polygon_2;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost            Cost;

int main()
{
  Polygon_2 polygon;
  std::cin >> polygon;
  Cost cost;
  polygon = PS::simplify(polygon, cost, Stop(0.5));
  
  std::cout.precision(12);
  std::cout << polygon << std::endl;

  return 0;
}


