#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/MP_Float.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitters.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>


typedef CGAL::Homogeneous<CGAL::MP_Float> R;

typedef R::Point_3 Point;

typedef R::FT FT;
typedef R::RT RT;

typedef CGAL::Plane_separator<FT> Separator;  
typedef CGAL::Kd_tree_traits_point<Point> Traits;
typedef CGAL::Creator_uniform_3<RT,Point> Creator; 

int main() {

  int bucket_size=10;
  
  const int data_point_number=1000; 
  typedef std::list<Point> point_list;
  point_list data_points;
  
 
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1000);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  
  Traits tr(bucket_size, 3, true);

  
  typedef CGAL::Kd_tree<Traits> Tree;


  Tree d(data_points.begin(), data_points.end(), tr);

  std::cout << "created kd tree using "  
  << data_point_number << " points. " << std::endl;
  d.statistics();

};


