// Approximate spatial searching: Example06.C
// Example illustrating for each separate splitting rule
// building a kd-tree   

#include <CGAL/basic.h>


#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/MP_Float.h>
#include <CGAL/Iso_rectangle_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/General_priority_search.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/L1_distance_rectangle_point.h>

typedef CGAL::Homogeneous<CGAL::MP_Float> R;
typedef CGAL::Point_3<R> Point;
typedef Point::R::FT FT;
typedef Point::R::RT RT;

typedef CGAL::Iso_cuboid_3<R> Rectangle;
typedef CGAL::Plane_separator<FT> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::L1_distance_rectangle_point<Rectangle,Point> L1_distance;
typedef CGAL::General_priority_search<Traits, Rectangle, L1_distance> 
NN_priority_search;

typedef CGAL::Creator_uniform_3<RT,Point> Creator;

int test_range_searching(CGAL::Split_rule_enumeration::Split_rule s) {

  int bucket_size=1;
  const int dim=3;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g(1000.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  
  Traits tr(bucket_size, s, 3.0, false);
  L1_distance tr_dist(dim);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // define range query
  // define range query
  int p[dim];
  int q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0;
        q[i]=  200;
  }

  Point P(p[0],p[1],p[2],1);  
  Point Q(q[0],q[1],q[2],1);
  Rectangle query_rectangle(P,Q);

 
  std::vector<NN_priority_search::Item_with_distance> nearest_neighbours;
  nearest_neighbours.reserve(nearest_neighbour_number);

  NN_priority_search NN(d, query_rectangle, tr_dist, 0);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  
  // NN.the_k_neighbours(std::back_inserter(nearest_neighbours));

  std::vector<NN_priority_search::Item_with_distance>::iterator it = nearest_neighbours.begin();

  CGAL::copy_n(NN.begin(), nearest_neighbour_number, it);
 
  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << *(nearest_neighbours[j].first) << std::endl; 
  }

  NN.statistics();
  
  return 0;
}; 
  
 

int main() {
  
  test_range_searching(CGAL::Split_rule_enumeration::MEDIAN_OF_MAX_SPREAD); 
  test_range_searching(CGAL::Split_rule_enumeration::MEDIAN_OF_RECTANGLE); 
  test_range_searching(CGAL::Split_rule_enumeration::MIDPOINT_OF_MAX_SPREAD);
  test_range_searching(CGAL::Split_rule_enumeration::MIDPOINT_OF_RECTANGLE);
  test_range_searching(CGAL::Split_rule_enumeration::FAIR);
  test_range_searching(CGAL::Split_rule_enumeration::SLIDING_MIDPOINT); 
  test_range_searching(CGAL::Split_rule_enumeration::SLIDING_FAIR);    

  return 0;
};


