#include <CGAL/basic.h>
#include <vector>
#include <iostream>

#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>

#include "../../examples/Spatial_searching/Point.h" 


typedef CGAL::Kd_tree_traits_point<double, Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef Tree::Splitter Splitter;
typedef CGAL::Creator_uniform_3<double,Point> Creator;

int main() {

  int bucket_size=10;
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));

  Splitter split(bucket_size);
  Tree d(data_points.begin(), data_points.end(), split);

  std::cout << "created kd tree using splitting rule "  
  << data_point_number << " points. " << std::endl;
  d.statistics(std::cout);
  return 0;

};


