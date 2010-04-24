// file: K_neighbor_search_with_circle.cpp
#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Euclidean_distance_sphere_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/point_generators_2.h>
#include <vector>
#include <algorithm>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef K::Circle_2 Circle;
typedef K::FT NT;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Euclidean_distance_sphere_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
  
int main() {

  const unsigned int N = 1000;
  const unsigned int K = 50;
  Random_points_iterator g(1.0);
  std::vector<Point> points;
  for(unsigned int i=0; i < N; i++){
    points.push_back(*g++);
  }
  
 
  Tree tree(points.begin(), points.end());

  double squared_radius = 0.04;

  Point center(0.2, 0.2);
  Circle query(center,squared_radius);

  Distance dist;

  Neighbor_search search(tree, query, K);
  std::vector<Point> result;

  double max_dist = -1;
  for (Neighbor_search::iterator it = search.begin();
       it != search.end();
       it++){ 
    result.push_back(it->first);
    if(CGAL::squared_distance(center, it->first) > max_dist){
      max_dist = it->second;
    }
  }
  assert(result.size() == K);
  std::sort(points.begin(), points.end());
  std::sort(result.begin(), result.end());
  
  std::vector<Point> diff;
  std::set_difference(points.begin(), points.end(),
		      result.begin(), result.end(),
		      std::back_inserter(diff));

  assert(diff.size() == N-K);
  for(std::vector<Point>::iterator it = diff.begin();
      it != diff.end();
      it++){
    assert(CGAL::squared_distance(center, *it) >= max_dist);
  }
  
  std::cout << "done" << std::endl;
  return 0;
} 
  
 


