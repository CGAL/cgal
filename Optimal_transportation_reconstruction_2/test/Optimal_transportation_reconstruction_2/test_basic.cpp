// test_basic.cpp

//----------------------------------------------------------
// Test the cgal environment for Optimal_transportation_reconstruction_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>
#include "testing_tools.h"

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                               FT;

int main ()
{
  std::vector<Point> points;
  //use the stair example for testing
  load_xy_file_points<Point>("data/stair-noise00.xy", points);

  CGAL::Optimal_transportation_reconstruction_2<K> otr2(points);

  if (otr2.run(100)) //100 steps
    std::cerr << "All done." << std::endl;
  else
    std::cerr << "Premature ending." << std::endl;

  otr2.print_stats_debug();
}
