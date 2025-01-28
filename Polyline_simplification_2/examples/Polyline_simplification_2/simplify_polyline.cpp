#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/IO/WKT.h>
#include <deque>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::deque<Point_2> Polyline_2;
typedef PS::Stop_above_cost_threshold Stop;
typedef PS::Squared_distance_cost Cost;
int main(int argc, char* argv[])
{
  Polyline_2 polyline;
  std::ifstream ifs( (argc==1)?"data/polyline.wkt":argv[1]);
  CGAL::IO::read_linestring_WKT(ifs, polyline);
  Cost cost;
  std::deque<Point_2> result;
  PS::simplify(polyline.begin(), polyline.end(), cost, Stop(0.5), std::back_inserter(result));

  std::cout.precision(12);
  for(std::size_t i=0; i < result.size(); ++i){
    std::cout << result[i] << std::endl;
  }
  return 0;
}
