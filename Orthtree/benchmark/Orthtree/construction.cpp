
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
  file.open((argc > 1) ? argv[1] : "../construction_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Build Time (ms) \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 10; num_points < 10000000; num_points *= 1.1) {

    // Create a collection of the right number of points
    auto points = generate<Kernel>(num_points);

    // Start the timer
    auto start = high_resolution_clock::now();

    // Build the tree
    Octree octree(points, points.point_map());
    octree.refine();

    // End the timer
    auto end = high_resolution_clock::now();

    file << num_points << ",";
    file << duration_cast<microseconds>(end - start).count() << "\n";

    std::cout << num_points << std::endl;
  }

  file.close();

  return 0;
}