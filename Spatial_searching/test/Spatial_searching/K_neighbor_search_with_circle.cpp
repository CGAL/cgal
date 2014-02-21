// file: K_neighbor_search_with_circle.cpp
#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Euclidean_distance_sphere_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Search_traits_adapter.h>
#include "Point_with_info.h"
#include <vector>
#include <algorithm>

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Circle_2 Circle;
typedef Kernel::FT NT;
typedef CGAL::Search_traits_2<Kernel> TreeTraits;
typedef CGAL::Euclidean_distance_sphere_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
//typdefs for Point_with_info
typedef Point_with_info_helper<Point>::type                                             Point_with_info;
typedef Point_property_map<Point>                                                       Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,TreeTraits>                   Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Distance>                         Distance_adapter;
typedef CGAL::K_neighbor_search<Traits_with_info, Distance_adapter>                   Neighbor_search_with_info;


const unsigned int N = 1000;
const unsigned int K = 50;

template <class K_search>
void run(std::vector<Point> points)
{
  typename K_search::Tree tree(
    boost::make_transform_iterator(points.begin(),Create_point_with_info<typename K_search::Point_d>()),
    boost::make_transform_iterator(points.end(),Create_point_with_info<typename K_search::Point_d>())
  );

  double squared_radius = 0.04;

  Point center(0.2, 0.2);
  Circle query(center,squared_radius);

  Distance dist;

  K_search search(tree, query, K);
  std::vector<Point> result;

  double max_dist = -1;
  for (typename K_search::iterator it = search.begin();
       it != search.end();
       it++){ 
    result.push_back(get_point(it->first));
    if(CGAL::squared_distance(center, get_point(it->first)) > max_dist){
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
}


int main() {


  Random_points_iterator g(1.0);
  std::vector<Point> points;
  for(unsigned int i=0; i < N; i++){
    points.push_back(*g++);
  }
  
  run<Neighbor_search>(points);
  run<Neighbor_search_with_info>(points);
 

  return 0;
} 
  
 


