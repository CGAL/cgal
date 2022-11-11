
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <chrono>
#include <cassert>

using namespace std::chrono;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>
Octree;

typedef CGAL::Search_traits_3<Kernel> Kd_tree_traits;
typedef CGAL::Orthogonal_k_neighbor_search<Kd_tree_traits> Kd_tree_search;
typedef Kd_tree_search::Tree Kd_tree;


void naive_vs_octree(std::size_t dataset_size) {

  std::cout << "[ " << dataset_size << " points ]" << std::endl;

  // Create a dataset
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(dataset_size);
  for (std::size_t i = 0; i < dataset_size; ++i)
    points.insert(*(generator++));

  // Choose another random point from the same bounds as the dataset
  Point random_point = *(generator++);

  // Use the naive algorithm to find the nearest point in the dataset
  Point naive_nearest = *points.points().begin();
  auto naive_start_time = high_resolution_clock::now();
  {

    FT distance_nearest = (std::numeric_limits<FT>::max)();
    for (auto &p : points.points()) {

      FT distance_current = CGAL::squared_distance(p, random_point);
      if (distance_current < distance_nearest) {

        distance_nearest = distance_current;
        naive_nearest = p;
      }
    }
  }
  duration<float> naive_elapsed_time = high_resolution_clock::now() - naive_start_time;

  std::cout << "Naive --> "
            << "distance^2 of "
            << CGAL::squared_distance(naive_nearest, random_point) << " "
            << "with a time of "
            << naive_elapsed_time.count() << "s "
            << std::endl;

  // Do the same using the octree
  Point octree_nearest = *generator;
  Octree octree(points, points.point_map());
  octree.refine(10, 20);
  auto octree_start_time = high_resolution_clock::now();
  {
    // TODO: Write a nearest-neighbor implementation and use it here
    std::vector<Point> k_neighbors;
    octree.nearest_neighbors(random_point, 1, std::back_inserter(k_neighbors));
    octree_nearest = *k_neighbors.begin();
  }
  duration<float> octree_elapsed_time = high_resolution_clock::now() - octree_start_time;

  std::cout << "Octree --> "
            << "distance^2 of "
            << CGAL::squared_distance(octree_nearest, random_point) << " "
            << "with a time of "
            << octree_elapsed_time.count() << "s "
            << std::endl;

  // Check that they produce the same answer
  assert(octree_nearest == naive_nearest);
}

void kdtree_vs_octree(std::size_t dataset_size, std::size_t K) {

  std::cout << "[ " << dataset_size << " points ]" << std::endl;

  // Create a dataset
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(dataset_size);
  for (std::size_t i = 0; i < dataset_size; ++i)
    points.insert(*(generator++));

  // Choose another random point from the same bounds as the dataset
  Point random_point = *(generator++);

  // Use the naive algorithm to find the nearest point in the dataset
  std::vector<Point> kd_tree_nearest_neighbors;
  Kd_tree kd_tree(points.points().begin(), points.points().end());
  kd_tree.build();
  auto kd_tree_start_time = high_resolution_clock::now();
  Kd_tree_search search(kd_tree, random_point, (unsigned int)(K));
  duration<float> kd_tree_elapsed_time = high_resolution_clock::now() - kd_tree_start_time;
  for (auto p : search)
    kd_tree_nearest_neighbors.push_back(p.first);

  std::cout << "Kd_tree --> "
            << kd_tree_nearest_neighbors.size() << " points "
            << "in " << kd_tree_elapsed_time.count() << "s "
            << std::endl;

  // Do the same using the octree
  std::vector<Point> octree_nearest_neighbors;
  Octree octree(points, points.point_map());
  octree.refine(10, 20);
  auto octree_start_time = high_resolution_clock::now();
  octree.nearest_neighbors(random_point, K, std::back_inserter(octree_nearest_neighbors));
  duration<float> octree_elapsed_time = high_resolution_clock::now() - octree_start_time;

  std::cout << "Octree --> "
            << octree_nearest_neighbors.size() << " points "
            << "in " << octree_elapsed_time.count() << "s "
            << std::endl;

  // Check that the octree produces the right number of results
  assert(octree_nearest_neighbors.size() == K);

  // Check that they produce the same answer
  for (std::size_t j = 0; j < K; ++j)
    assert(octree_nearest_neighbors[j] == kd_tree_nearest_neighbors[j]);

}

int main(void) {

  naive_vs_octree(500);
  naive_vs_octree(1000);
  naive_vs_octree(10000);
  naive_vs_octree(100000);

  kdtree_vs_octree(100, 16);
  kdtree_vs_octree(1000, 16);
  kdtree_vs_octree(10000, 16);
  kdtree_vs_octree(100000, 16);

  return EXIT_SUCCESS;
}
