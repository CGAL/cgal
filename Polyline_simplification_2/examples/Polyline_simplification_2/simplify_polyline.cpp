#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <list>
#include <vector>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::list<Point_2> Polyline_2;
typedef PS::Stop_above_cost_threshold Stop;
typedef PS::Squared_distance_cost Cost;

int main()
{
  Polyline_2 polyline;
  int n;
  Point_2 p;
  std::cin >> n;
  while(std::cin >> p){
    polyline.push_back(p);
  }
  Cost cost;
  std::vector<Point_2> result;
  PS::simplify(polyline.begin(), polyline.end(), cost, Stop(0.5), std::back_inserter(result));
  
  std::cout.precision(12);
  for(std::size_t i=0; i < result.size(); ++i){
    std::cout << result[i] << std::endl;
  }
  return 0;
}


