#include <boost/config.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/IO/WKT.h>


namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                   Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>        Polygon_with_holes_2;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost            Cost;

int main(int argc, char* argv[])
{
  std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);
  Polygon_with_holes_2 polygon;
  CGAL::read_polygon_WKT(ifs, polygon);
  Cost cost;
  polygon = PS::simplify(polygon, cost, Stop(0.5));

  std::cout.precision(12);
  CGAL::write_polygon_WKT(std::cout, polygon) << std::endl;

  return 0;
}
#else
int main()
{
  return 0;
}
#endif
