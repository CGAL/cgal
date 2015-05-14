//! \file examples/Minkowski_sum_2/sum_by_decomposition.cpp
// Computing the Minkowski sum of two non-convex polygons read from a file
// using the small-side angle-bisector decomposition strategy.

#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>

#include "pgn_print.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;

int main(int argc, char* argv[])
{
  // Open the input file and read two polygons from it.
  const char* filename = (argc > 1) ? argv[1] : "rooms_star.dat";
  std::ifstream    in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }
  Polygon_2   P, Q;
  in_file >> P >> Q;
  in_file.close();

  // Compute the Minkowski sum using the decomposition approach.
  CGAL::Small_side_angle_bisector_decomposition_2<Kernel>  ssab_decomp;
  Polygon_with_holes_2  sum = CGAL::minkowski_sum_2(P, Q, ssab_decomp);
  std::cout << "P (+) Q = "; print_polygon_with_holes(sum);
  return 0;
}
