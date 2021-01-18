
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

  int num_runs = 100;

  size_t k = 10;

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "../nearest_neighbor_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Octree,kDTree \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 100; num_points < 100000; num_points *= 1.05) {

    // We want the average of several runs for each point count, for cleaner results
    float octreeAverage = 0;
    float kdtreeAverage = 0;
    float naiveAverage = 0;

    // Repeat the tests, generating a new point set for each run
    for (int i = 0; i < num_runs; ++i) {

      // Create a collection of the right number of points
      auto points = generate<Kernel>(num_points);

      // Create a search point
      auto search_point = *(generate<Kernel>().points().end() - 1);

      // Build the kd tree from the point set
      Kdtree kdtree(points.points().begin(), points.points().end());
      kdtree.build();

      // Time how long it takes to find neighbors using the kd tree
      auto kdtreeTime = bench<microseconds>(
              [&] {
                Kd_tree_search search(kdtree, search_point, k);
              }
      );

      // Time how long it takes to find neighbors using a naive approach
//    auto naiveTime = bench<microseconds>(
//            [&] {
//
//              std::vector<Point> nearest_neighbors;
//
//              // Iterate over every point
//              for (auto &point : points.points()) {
//
//                // Find out how this point ranks in comparison with other points we've saved
//                auto iter = nearest_neighbors.begin();
//                for (; iter < nearest_neighbors.end() &&
//                       CGAL::squared_distance(point, search_point) <
//                       CGAL::squared_distance(*iter, search_point);
//                       iter++) {}
//
//                // Add the point to the list (it'll usually be at the end)
//                nearest_neighbors.insert(iter, point);
//
//                // Never keep more than k neighbors
//                if (nearest_neighbors.size() > k)
//                  nearest_neighbors.resize(k);
//
//              }
//
//            }
//    );

      // Build the octree from points (this had to be done second because it rearranges the point set)
      Octree octree(points, points.point_map());
      octree.refine();

      // Time how long it takes to find neighbors using the octree
      auto octreeTime = bench<microseconds>(
              [&] {
                std::vector<Point> nearest_neighbors;
                octree.nearest_neighbors(search_point, k, std::back_inserter(nearest_neighbors));
              }
      );

      // Incorporate our results into the average
      octreeAverage += (float) octreeTime.count() / (float) num_runs;
      kdtreeAverage += (float) kdtreeTime.count() / (float) num_runs;
//      naiveAverage += (float) naiveTime.count() / (float) num_runs;

      // A simple progress indication
      std::cout << ".";
    }

    file << num_points << ",";
    file << octreeAverage << ",";
    file << kdtreeAverage << ",";
//    file << naiveAverage << ",";
    file << "\n";

    std::cout << num_points << std::endl;

  }

  file.close();

  return 0;
}
