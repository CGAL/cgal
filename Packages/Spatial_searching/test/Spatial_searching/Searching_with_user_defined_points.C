#include <CGAL/basic.h>

#include <vector>
#include <cassert>
#include <iostream>

#include <CGAL/Kd_tree.h>
#include <CGAL/Random.h>
#include <CGAL/Splitters.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>

#ifdef CUSTOM_POINTS

#include <CGAL/Kd_tree_traits_point.h>
#include "../../examples/Spatial_searching/Point.h" 
typedef CGAL::Kd_tree_traits_point<double, Point, const double*, Construct_coord_iterator> Traits;

#else

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree_traits_point_3.h>
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Kd_tree_traits_point_3<K> Traits;
typedef CGAL::Euclidean_distance<Traits> Distance;
#endif

typedef CGAL::Creator_uniform_3<double,Point> Creator;

typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> NN_orthogonal_search;
typedef CGAL::K_neighbor_search<Traits, Distance> NN_general_search;
typedef NN_general_search::Tree Tree;
typedef NN_general_search::Splitter Splitter;

typedef std::vector<Traits::Point> Vector;
typedef std::vector<Point> Query_vector;

int main() {

  int bucket_size=10;
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  
  Splitter split1(bucket_size);

  NN_orthogonal_search::Tree d1(data_points.begin(), data_points.end(), split1);

  
  
  std::cout << "created kd tree using extended nodes containing "   
  << data_point_number << " points. " << std::endl;
  d1.statistics(std::cout);

  
  Splitter split2(bucket_size);
  NN_general_search::Tree d2(data_points.begin(), data_points.end(), split2);

  std::cout << "created kd tree using no extended nodes containing "
  << data_point_number << " points. " << std::endl;
  d2.statistics(std::cout);
  
  // neighbour searching
  const int query_point_number=5;
  Query_vector query_points;
  CGAL::copy_n( g, query_point_number+1, std::back_inserter(query_points));
  
  // nearest neighbour searching using extended nodes
  std::vector<NN_orthogonal_search::Point_with_distance> nearest_neighbours1;
  // nearest_neighbours1.reserve(query_point_number+1);
  
  
  // nearest neighbour searching using no extended nodes
  std::vector<NN_general_search::Point_with_distance> nearest_neighbours2;
  // nearest_neighbours2.reserve(query_point_number+1);
  
  for (int i=1; i < query_point_number+1; ++i) { 
     NN_orthogonal_search NN1(d1, query_points[i], 1, 0.0);
     std::cout << "neighbour searching statistics using extended nodes: " << std::endl;
     NN1.statistics(std::cout);
     NN1.the_k_neighbors(std::back_inserter(nearest_neighbours1));
     NN_general_search NN2(d2, query_points[i], 1, 0.0, false);
     std::cout << "neighbor searching statistics using no extended nodes: " << std::endl;
     NN2.statistics(std::cout);
     NN2.the_k_neighbors(std::back_inserter(nearest_neighbours2));
  }
  
  std::cout << "results neighbor searching:" << std::endl;

  for (int j=0; j < query_point_number; ++j) { 
     std::cout << " d(q,nearest neighbor)=" << nearest_neighbours1[j].second << 
     " d(q,furthest neighbour)=" << nearest_neighbours2[j].second << std::endl; 
  } 

  return 0;
};





