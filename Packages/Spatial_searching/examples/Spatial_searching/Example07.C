// Approximate spatial searching: Example07.C
// Example illustrating for each separate splitting rule
// building a kd-tree 
#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Random.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/L1_distance_rectangle_point.h>
#include <CGAL/Weighted_Minkowski_distance.h>

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::R::FT NT;

typedef CGAL::Kd_tree_rectangle<NT> Rectangle;
typedef CGAL::Plane_separator<NT> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Weighted_Minkowski_distance<Point, Point> Distance;
typedef CGAL::Orthogonal_standard_search<Traits, Point, Distance> 
NN_standard_search;
  

int test_range_searching(CGAL::Split_rule_enumeration::Split_rule s) {

  int bucket_size=1;
  const int dim=4;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) {
	    NT v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  
  Traits tr(bucket_size, s, 3.0, false);

  Distance::Weight_vector w(4);
  w[0]=1.0; w[1]=1.0; w[2]=1.0; w[3]=1.0;

  Distance tr_dist(2,dim,w);

  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // define query item
  double q[dim];
  for (int i=0; i<dim; i++) {
  	q[i]=0.5;
  }
  Point query_item(dim,q,q+dim);

  std::vector<NN_standard_search::Item_with_distance> nearest_neighbours;
  nearest_neighbours.reserve(nearest_neighbour_number);

  NN_standard_search NN(d, query_item, tr_dist, nearest_neighbour_number, 0.0);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  NN.statistics();
  NN.the_k_neighbours(std::back_inserter(nearest_neighbours));

  for (int i=0; i < nearest_neighbour_number; ++i) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[i].second << 
     " nn= " << *(nearest_neighbours[i].first) << std::endl; 
  }
  
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


