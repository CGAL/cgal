#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/MP_Float.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitters.h>
#include <CGAL/General_priority_search.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/Manhattan_distance_rectangle_point.h>

typedef CGAL::Homogeneous<CGAL::MP_Float> R;

typedef CGAL::Point_3<R> Point;
typedef Point::R::FT FT;
typedef Point::R::RT RT;

typedef CGAL::Iso_cuboid_3<R> Rectangle;
typedef CGAL::Plane_separator<FT> Separator;

typedef CGAL::Kd_tree_traits_point<Point> Traits;
typedef CGAL::Manhattan_distance_rectangle_point<Rectangle,Point> L1_distance;
typedef CGAL::General_priority_search<Traits, L1_distance, Rectangle> 
NN_priority_search;

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
  
  Traits tr(bucket_size, 3.0, false);
  L1_distance tr_dist(dim);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

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

  NN_priority_search NN(d, query_rectangle, tr_dist, 0);
  std::cout << "neighbour searching statistics without using extended nodes: " << std::endl;
  
  
  int neighbour_number_copy=nearest_neighbour_number;
  
  
  CGAL::copy_n(NN.begin(), neighbour_number_copy, std::back_inserter(nearest_neighbours));
  
 
  
  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << *(nearest_neighbours[j].first) << std::endl; 
  }

  NN.statistics();
  
  return 0;
}; 



