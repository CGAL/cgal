#include <CGAL/Homogeneous.h>
#include <vector>
#include <cassert>

#include <iostream>

#include <CGAL/MP_Float.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitters.h>
#include <CGAL/Incremental_neighbor_search.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>

typedef CGAL::Homogeneous<CGAL::MP_Float> R;

typedef R::Point_3 Point;
typedef R::FT FT;
typedef R::RT RT;

typedef R::Iso_cuboid_3 Rectangle;

typedef CGAL::Kd_tree_traits_point_3<R> Traits;
typedef CGAL::Manhattan_distance_iso_box_point<Traits, Rectangle> L1_distance;
typedef CGAL::Incremental_neighbor_search<Traits, L1_distance>  NN_priority_search;
typedef NN_priority_search::Tree Tree;
typedef NN_priority_search::Splitter Splitter;

typedef CGAL::Creator_uniform_3<RT,Point> Creator;

int main() {
  
  int bucket_size=1;
  const int dim=3;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g(1000.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
Splitter split(bucket_size);

  Tree d(data_points.begin(), data_points.end(), split);

  double p[dim];
  double q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0.0;
        q[i]=  600.0;
  }

  Point P(p[0],p[1],p[2]);  
  Point Q(q[0],q[1],q[2]);

  Rectangle query_rectangle(P,Q);

  
  std::vector<NN_priority_search::Point_with_distance> nearest_neighbours;

  NN_priority_search NN(d, query_rectangle, 0);
  std::cout << "neighbour searching statistics without using extended nodes: " << std::endl;
  
  
  int neighbour_number_copy=nearest_neighbour_number;
  
  
  CGAL::copy_n(NN.begin(), neighbour_number_copy, std::back_inserter(nearest_neighbours));
  
 
  
  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << nearest_neighbours[j].first << std::endl; 
  }

  NN.statistics(std::cout);

  return 0;
} 



