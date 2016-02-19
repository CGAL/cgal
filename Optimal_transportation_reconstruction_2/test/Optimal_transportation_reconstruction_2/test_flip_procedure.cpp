
// test_flip_procedure.cpp

//----------------------------------------------------------
// Test the cgal environment for Optimal_transportation_reconstruction_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>
#include "testing_tools.h"

#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                               FT;

int main ()
{
  std::vector<Point> points;
  //use the stair example for testing
  load_xy_file_points<Point>("data/stair-noise00.xy", points);

  for (std::size_t i = 1 ; i <= points.size() ; i += 20)
  {
    CGAL::Optimal_transportation_reconstruction_2<K> otr2(points);
    otr2.run_until(i);
    otr2.print_stats_debug();
    assert(otr2.number_of_vertices() == i);
  }
}
