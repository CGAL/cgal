//! \file examples/Arrangement_on_surface_2/predefined_kernel.cpp
// Constructing an arrangement of intersecting line segments using the
// predefined kernel with exact constructions and exact predicates.

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include "arr_exact_construction_segments.h"
#include "arr_print.h"
#include "read_objects.h"

int main (int argc, char* argv[]) {
  // Get the name of the input file from the command line, or use the default
  // fan_grids.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "fan_grids.dat";

  // Open the input file.
  std::ifstream in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << " ...\n";
    return 1;
  }

  std::list<Segment> segments;
  read_objects<Segment>(filename, std::back_inserter(segments));

  // Construct the arrangement by aggregately inserting all segments.
  Arrangement arr;
  CGAL::Timer timer;

  std::cout << "Performing aggregated insertion of "
            << segments.size() << " segments.\n";

  timer.start();
  insert(arr, segments.begin(), segments.end());
  timer.stop();

  print_arrangement_size(arr);

  std::cout << "Construction took " << timer.time() << " seconds.\n";

  return 0;
}
