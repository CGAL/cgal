// Approximate spatial searching: Example05.C
// Example illustrating for each separate splitting rule
// building a kd-tree 

#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/Iso_rectangle_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Random.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/General_standard_search.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/L1_distance_rectangle_point.h>

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::R::FT NT;

typedef CGAL::Iso_rectangle_d<R> Rectangle;
typedef CGAL::Plane_separator<NT> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::L1_distance_rectangle_point<Rectangle,Point> L1_distance;
typedef CGAL::General_standard_search<Traits, Rectangle, L1_distance> 
NN_standard_search;
  

int test_range_searching(CGAL::Split_rules::Split_rule s) {

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
  L1_distance tr_dist(dim);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // define range query
  int p[dim];
  int q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0;
        q[i]=  1000;
  }

  Point P(p[0],p[1],p[2],1000);
  Point Q(q[0],q[1],q[2],1000);

  Rectangle query_rectangle(P,Q);

  std::vector<NN_standard_search::Item_with_distance> nearest_neighbours;
  
  NN_standard_search NN(d, query_rectangle, tr_dist, nearest_neighbour_number, 0.0);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  NN.statistics();
  NN.the_k_neighbours(std::back_inserter(nearest_neighbours));

  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << *(nearest_neighbours[j].first) << std::endl; 
  }
  
  return 0;
}; 
  
 

int main() {
  
  test_range_searching(CGAL::Split_rules::MEDIAN_OF_MAX_SPREAD); 
  test_range_searching(CGAL::Split_rules::MEDIAN_OF_RECTANGLE); 
  test_range_searching(CGAL::Split_rules::MIDPOINT_OF_MAX_SPREAD);
  test_range_searching(CGAL::Split_rules::MIDPOINT_OF_RECTANGLE);
  test_range_searching(CGAL::Split_rules::FAIR);
  test_range_searching(CGAL::Split_rules::SLIDING_MIDPOINT); 
  test_range_searching(CGAL::Split_rules::SLIDING_FAIR);    

  return 0;
};


