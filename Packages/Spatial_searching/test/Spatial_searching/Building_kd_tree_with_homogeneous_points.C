
#include <CGAL/Homogeneous.h>

#include <CGAL/MP_Float.h>
#include <iostream>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>


typedef CGAL::Homogeneous<CGAL::MP_Float> R;

typedef R::Point_3 Point;

typedef R::FT FT;
typedef R::RT RT;


typedef CGAL::Kd_tree_traits_point_3<R> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef Tree::Splitter Splitter;
typedef CGAL::Creator_uniform_3<RT,Point> Creator; 

int main() {

  int bucket_size=10;
  
  const int data_point_number=1000; 
  typedef std::list<Point> point_list;
  point_list data_points;
  
 
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1000);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  
  Splitter  split(bucket_size);
  
  Tree d(data_points.begin(), data_points.end(), split);

  std::cout << "created kd tree using "  
  << data_point_number << " points. " << std::endl;
  d.statistics(std::cout);

};


