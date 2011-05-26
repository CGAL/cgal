// file: Building_kd_tree_with_ddim_points.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/point_generators_d.h>
#include <vector>
#include <cassert>
#include <iostream>

//typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Homogeneous_d<double> K;
typedef CGAL::Point_d<K> Point;

typedef CGAL::Random_points_in_cube_d<Point> Point_generator;

typedef K SearchTraits;
typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits> Orthogonal_k_neighbor_search;
typedef CGAL::K_neighbor_search<SearchTraits>  K_neighbor_search;


int main() {
  
  const unsigned int N = 1000;
  const unsigned int K = 20;
  
  std::vector<Point> points;
  
  Point_generator g(3);
  CGAL::cpp0x::copy_n( g, N, std::back_inserter(points));
  g++;
  Point query = *g;
  
  Orthogonal_k_neighbor_search::Tree o_tree(points.begin(), points.end());
  o_tree.statistics(std::cout);

  K_neighbor_search::Tree tree(points.begin(), points.end());
  tree.statistics(std::cout);

  // do checking
  double dist = 0;
  std::vector<Point> result;
  Orthogonal_k_neighbor_search o_search(o_tree, query,K);
  for(Orthogonal_k_neighbor_search::iterator it = o_search.begin();
      it != o_search.end();
      it++){
    result.push_back(it->first);
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      assert(CGAL::squared_distance(query, *it) >= dist);
    }
  }

  result.clear();
  K_neighbor_search search(tree, query,K);
  for(K_neighbor_search::iterator it = search.begin();
      it != search.end();
      it++){
    result.push_back(it->first);
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      assert(CGAL::squared_distance(query, *it) >= dist);
    }
  }
  std::cout << "done" << std::endl;
  return 0;
}


