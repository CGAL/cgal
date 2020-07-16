
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>

#include <chrono>
#include <cassert>

using namespace std::chrono;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
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

    FT distance_nearest = std::numeric_limits<FT>::max();
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
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << naive_nearest << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(naive_nearest, random_point) << " "
            << "with a time of "
            << naive_elapsed_time.count()
            << std::endl;

  // Do the same using the octree
  Point octree_nearest = *generator;
  auto point_map = points.point_map();
  Octree octree(points, point_map);
  octree.refine(10, 1);
  auto octree_start_time = high_resolution_clock::now();
  {
    // TODO: Write a nearest-neighbor implementation and use it here
    std::vector<Point> k_neighbors;
    octree.nearest_k_neighbours(random_point, 1, std::back_inserter(k_neighbors));
    octree_nearest = *k_neighbors.begin();
  }
  duration<float> octree_elapsed_time = high_resolution_clock::now() - octree_start_time;

  std::cout << "Octree --> "
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << octree_nearest << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(octree_nearest, random_point) << " "
            << "with a time of "
            << octree_elapsed_time.count()
            << std::endl;

  // Check that they produce the same answer
  assert(octree_nearest == naive_nearest);

  // Check that the octree was faster
  assert(octree_elapsed_time < naive_elapsed_time);
}

void kdtree_vs_octree(std::size_t dataset_size) {

  std::cout << "[ " << dataset_size << " points ]" << std::endl;

  int K = 16;

  // Create a dataset
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(dataset_size);
  for (std::size_t i = 0; i < dataset_size; ++i)
    points.insert(*(generator++));

  // Choose another random point from the same bounds as the dataset
  Point random_point = *(generator++);

  // Use the naive algorithm to find the nearest point in the dataset
  std::vector<Point> kd_tree_nearest_neighbours;
  Kd_tree kd_tree(points.points().begin(), points.points().end());
  {
    Kd_tree_search search(kd_tree, random_point, K);
    for (auto p : search)
      kd_tree_nearest_neighbours.push_back(p.first);
  }

  std::cout << "Kd_tree --> "
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << kd_tree_nearest_neighbours[0] << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(kd_tree_nearest_neighbours[0], random_point)
            << std::endl;

  // Do the same using the octree
  std::vector<Point> octree_nearest_neighbours;
  auto point_map = points.point_map();
  Octree octree(points, point_map);
  octree.refine(10, 1);
  {
    octree.nearest_k_neighbours(random_point, K, std::back_inserter(octree_nearest_neighbours));
  }

  std::cout << "Octree --> "
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << octree_nearest_neighbours[0] << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(octree_nearest_neighbours[0], random_point)
            << std::endl;

  // Check that they produce the same answer
  for (int j = 0; j < K; ++j) {

    std::cout << j << std::endl;
    assert(octree_nearest_neighbours[j] == kd_tree_nearest_neighbours[j]);
  }
}

int main(void) {

//  naive_vs_octree(100);
//  naive_vs_octree(1000);
  naive_vs_octree(10000);
  naive_vs_octree(100000);

  kdtree_vs_octree(100);
//  kdtree_vs_octree(1000);
//  kdtree_vs_octree(10000);
//  kdtree_vs_octree(100000);

  return EXIT_SUCCESS;
}
