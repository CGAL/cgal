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

// #include <CGAL/MP_Float.h>
#include <CGAL/Iso_rectangle_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/General_priority_search.h>
// #include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/L1_distance_rectangle_point.h>

// typedef CGAL::Homogeneous<CGAL::MP_Float> R;
typedef CGAL::Cartesian<double> R;
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

template <class InputIterator, class Size, class OutputIterator>
OutputIterator my_copy_n( InputIterator first, Size n,
                       OutputIterator result) {
  
  Size number_of_el_to_compute=n;
  
  while( (number_of_el_to_compute > 0)) {
    
    number_of_el_to_compute--;
    
    *result = *first;
    first++;
    result++;
  }
 
  return result;
}


int test_range_searching(CGAL::Split_rules::Split_rule s) {

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
  /*
  int p[dim];
  int q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0;
        q[i]=  200;
  }

  Point P(p[0],p[1],p[2],1);  
  Point Q(q[0],q[1],q[2],1);
  */

  
  double p[dim];
  double q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0.0;
        q[i]=  200.0;
  }

  Point P(p[0],p[1],p[2]);  
  Point Q(q[0],q[1],q[2]);

  Rectangle query_rectangle(P,Q);

 
  std::vector<NN_priority_search::Item_with_distance> nearest_neighbours;
  nearest_neighbours.reserve(nearest_neighbour_number);

  NN_priority_search NN(d, query_rectangle, tr_dist, 0);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  
  

  std::vector<NN_priority_search::Item_with_distance>::iterator it = nearest_neighbours.begin();

  // CGAL::copy_n(NN.begin(), nearest_neighbour_number, it);
  // NN_priority_search::iterator NN_it=NN.begin();
  // use non const int neighbour_number instead of const int nearest_neighbour nubmer
  // to avoid Borland compiler bug
  int neighbour_number_copy=nearest_neighbour_number;
  CGAL::copy_n(NN.begin(), neighbour_number_copy, it);
  // my_copy_n(NN_it, neighbour_number, it);
  // my_copy_n(NN.begin(), neighbour_number, it);


  // my_copy_n(NN.begin(), nearest_neighbour_number, it); 
  // compiler  i686_CYGWINNT-5.0-1.3.2_bcc32.exe-5.51 complained
  // cannot modify constant object in function my_copy_n

  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << *(nearest_neighbours[j].first) << std::endl; 
  }

  NN.statistics();
  
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


