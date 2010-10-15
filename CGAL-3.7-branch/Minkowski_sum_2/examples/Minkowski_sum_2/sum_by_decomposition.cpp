//! \file examples/Minkowski_sum_2/sum_by_decomposition.cpp
// Computing the Minkowski sum of two non-convex polygons read from a file
// using the small-side angle-bisector decomposition strategy.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <iostream>
#include <fstream>

#include "print_utils.h"

struct Kernel : public CGAL::Exact_predicates_exact_constructions_kernel {};

typedef Kernel::Point_2                               Point_2;
typedef CGAL::Polygon_2<Kernel>                       Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>            Polygon_with_holes_2;

int main ()
{
  // Open the input file.
  std::ifstream    in_file ("rooms_star.dat");

  if (! in_file.is_open())
  {
    std::cerr << "Failed to open the input file." << std::endl;
    return (1);
  }

  // Read the two polygons from the file and compute their Minkowski sum.
  Polygon_2   P, Q;

  in_file >> P >> Q;
  in_file.close();

  // Compute the Minkowski sum using the decomposition approach.
  CGAL::Small_side_angle_bisector_decomposition_2<Kernel>  ssab_decomp;

  Polygon_with_holes_2  sum = minkowski_sum_2 (P, Q, ssab_decomp);

  std::cout << "P (+) Q = "; print_polygon_with_holes (sum);

  return (0);
}
