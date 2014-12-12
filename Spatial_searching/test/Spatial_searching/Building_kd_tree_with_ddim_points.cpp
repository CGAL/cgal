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
#include <CGAL/Search_traits_adapter.h>
#include "Point_with_info.h"
#include <vector>
#include <cassert>
#include <iostream>

//typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Homogeneous_d<double> Kernel;
typedef CGAL::Point_d<Kernel> Point;

typedef CGAL::Random_points_in_cube_d<Point> Point_generator;

typedef Kernel SearchTraits;
typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits>        Orthogonal_k_neighbor_search;
typedef CGAL::K_neighbor_search<SearchTraits>                   K_neighbor_search;
typedef Orthogonal_k_neighbor_search::Distance                  Distance;
//typdefs for Point_with_info
typedef Point_with_info_helper<SearchTraits::Point_d>::type                             Point_with_info;
typedef Point_property_map<SearchTraits::Point_d>                                     Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,SearchTraits>               Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Distance>                       Distance_adapter;
typedef CGAL::Orthogonal_k_neighbor_search<Traits_with_info,Distance_adapter>         Orthogonal_k_neighbor_search_with_info;
typedef CGAL::K_neighbor_search<Traits_with_info,Distance_adapter>                    K_neighbor_search_with_info;


const unsigned int N = 1000;
const unsigned int K = 20;

template <class OK_search,class K_search>
void run(const std::vector<Point>& points,const Point& query)
{
  typename OK_search::Tree o_tree(
      boost::make_transform_iterator(points.begin(),Create_point_with_info<typename OK_search::Point_d>()), 
      boost::make_transform_iterator(points.end(),Create_point_with_info<typename OK_search::Point_d>())
  );
  
  o_tree.statistics(std::cout);

  typename K_search::Tree tree(
      boost::make_transform_iterator(points.begin(),Create_point_with_info<typename K_search::Point_d>()), 
      boost::make_transform_iterator(points.end(),Create_point_with_info<typename K_search::Point_d>())
  );
  tree.statistics(std::cout);

  // do checking
  double dist = 0;
  std::vector<Point> result;
  OK_search o_search(o_tree, query,K);
  for(typename OK_search::iterator it = o_search.begin();
      it != o_search.end();
      it++){
    result.push_back(get_point(it->first));
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::const_iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      assert(CGAL::squared_distance(query, *it) >= dist);
    }
  }

  result.clear();
  K_search search(tree, query,K);
  for(typename K_search::iterator it = search.begin();
      it != search.end();
      it++){
    result.push_back(get_point(it->first));
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::const_iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      assert(CGAL::squared_distance(query, *it) >= dist);
    }
  }
  std::cout << "done" << std::endl;  
}

int main() {
  
 
  std::vector<Point> points;
  
  Point_generator g(3);
  CGAL::cpp11::copy_n( g, N, std::back_inserter(points));
  g++;
  Point query = *g;

  run<Orthogonal_k_neighbor_search,K_neighbor_search>(points,query);
  run<Orthogonal_k_neighbor_search_with_info,K_neighbor_search_with_info>(points,query);

  return 0;
}


