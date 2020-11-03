
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

int main(int argc, char **argv) {

  size_t k = 10;

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "../nearest_neighbor_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Octree,kDTree \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 100; num_points < 10000000; num_points *= 1.05) {

    file << num_points << ",";
    std::cout << num_points << std::endl;

    // Create a collection of the right number of points
    auto points = generate<Kernel>(num_points);

    // Create a search point
    auto search_point = *generate<Kernel>().points().begin();

    {
      // Build the tree (not timed)
      auto octreePoints = points;
      Octree octree(octreePoints, octreePoints.point_map());
      octree.refine();

      auto octreeTime = bench<microseconds>(
              [&] {
                std::vector<Point> nearest_neighbors;
                octree.nearest_neighbors(search_point, k, std::back_inserter(nearest_neighbors));
              }
      );
      file << octreeTime.count() << ",";
    }

    {
      auto kdtreePoints = points;
      Kdtree kdtree(kdtreePoints.points().begin(), kdtreePoints.points().end());
      kdtree.build();

      auto kdtreeTime = bench<microseconds>(
              [&] {
                Kd_tree_search search(kdtree, search_point, k);
              }
      );
      file << kdtreeTime.count() << "\n";
    }

  }

  file.close();

  return 0;
}
