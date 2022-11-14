
#define CGAL_TRACE_STREAM std::cerr

#include "util.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Octree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <iostream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;

typedef CGAL::Search_traits_3<Kernel> Kd_tree_traits;
typedef CGAL::Orthogonal_k_neighbor_search<Kd_tree_traits> Kd_tree_search;
typedef Kd_tree_search::Tree Kdtree;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::milliseconds;

int main(int argc, char **argv) {

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "../construction_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Octree,kDTree \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 10; num_points < 10000000; num_points *= 1.1) {

    // Create a collection of the right number of points
    auto points = generate<Kernel>(num_points);

    auto octreePoints = points;
    auto octreeTime = bench<milliseconds>(
            [&] {
              // Build the tree
              Octree octree(octreePoints, octreePoints.point_map());
              octree.refine();
            }
    );

    auto kdtreePoints = points;
    auto kdtreeTime = bench<milliseconds>(
            [&] {
              // Build the tree
              Kdtree kdtree(kdtreePoints.points().begin(), kdtreePoints.points().end());
              kdtree.build();
            }
    );

    file << num_points << ",";
    file << octreeTime.count() << ",";
    file << kdtreeTime.count() << "\n";

    std::cout << num_points << std::endl;
  }

  file.close();

  return 0;
}
