// test_reconstruction_until.cpp

//----------------------------------------------------------
// Test the cgal environment for Optimal_transportation_reconstruction_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                               FT;
typedef K::Segment_2                                         Segment;

int main ()
{
  std::vector<Point> points;

  //use the stair example for testing
  load_xy_file_points<Point>("data/stair-noise00.xy", points);

  CGAL::Optimal_transportation_reconstruction_2<K> otr2(points);
  otr2.run_until(9);
  otr2.print_stats_debug();

  std::vector<Point> isolated_points;
  std::vector<Segment> edges;

  otr2.list_output(
    std::back_inserter(isolated_points), std::back_inserter(edges));

  std::cout << "Isolated_points: " << isolated_points.size() << std::endl;
  std::cout << "Edges: " << edges.size() << std::endl;

  assert(isolated_points.size() == 0);
  assert(edges.size() == 8);
}
