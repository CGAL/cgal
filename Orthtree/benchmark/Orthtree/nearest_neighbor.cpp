
#define CGAL_TRACE_STREAM std::cerr

#include "util.h"

#include <CGAL/Octree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;

int main(int argc, char **argv) {

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "../nearest_neighbor_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Search Time (ms) \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 10; num_points < 1000000; num_points *= 1.1) {

    // Create a collection of the right number of points
    auto points = generate<Kernel>(num_points);

    // Create a search point
    auto search_point = *generate<Kernel>().points().begin();

    // Build the tree (not timed)
    Octree octree(points, points.point_map());
    octree.refine();

    // Start the timer
    auto start = high_resolution_clock::now();

    // Find the nearest point to the search point
    std::vector<Point> nearest_neighbors;
    octree.nearest_neighbors(search_point, 10, std::back_inserter(nearest_neighbors));

    // End the timer
    auto end = high_resolution_clock::now();

    file << num_points << ",";
    file << duration_cast<microseconds>(end - start).count() << "\n";
  }

  file.close();

  return 0;
}
