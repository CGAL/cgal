// example using nearest_neighbour_iterator for L2
// benchmark example using 10000 data points and 2000 query points
// both generated with Random_points_in_cube_3<Point_3>

#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>

#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_L2.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>

int test_benchmark_nearest_neighbour_L2() {

typedef CGAL::Cartesian<double> R;
typedef R::Point_3 Point;
typedef CGAL::Creator_uniform_3<double,Point> Creator;
typedef std::vector<Point> Vector;

typedef CGAL::Kernel_traits<Point>::Kernel K;
typedef K::FT NT;

typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Nearest_neighbour_L2<Traits,CGAL::Search_nearest_neighbour>::iterator NNN_Iterator;

  CGAL::Timer t;
  int dim=3;
  int point_number=10000;
  int query_point_number=2000;
  int bucket_size=1;
  double eps=0.0;

  std::cout << "test parameters: d=" << dim << " point_number=" << point_number << std::endl;
  std::cout << "query_point_number=" << query_point_number << " bucket_size="
  << bucket_size << " eps=" << eps << std::endl;

  // generate 10000 data points
  // Prepare a vector for 10000 points.
    
  // Vector data_points;
  // data_points.reserve(point_number);

  typedef std::list<Point> point_list;
  point_list data_points;
  
  // Create 10000 points within a cube.
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, point_number, std::back_inserter(data_points));

  // generate 10000 data points
  // Prepare a vector for 10000 points.
    
  
  Vector query_points;
  query_points.reserve(2000);

  // Create 2000 query points within the same cube.
  CGAL::copy_n( g, query_point_number, std::back_inserter(query_points));

  t.reset(); t.start();
  
  Traits tr(bucket_size, CGAL::SLIDING_MIDPOINT, 3.0, true);
  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr, true);
  t.stop();

  std::cout << "created binary search tree containing" << std::endl
  << point_number << " random points in the 3-dim unit cube in time " << t.time() <<
  " seconds " << std::endl;
  d.statistics();

  // end of building binary search tree


  /*
  std::vector<double> v(dim);
  */
  
  std::vector<Traits::Item_with_distance> nearest_neighbours(query_point_number);

  t.reset(); t.start();
  for (int i=0; i < query_point_number; i++) {
  
    NNN_Iterator NNN_Iterator1(d,query_points[i],0.0);
    nearest_neighbours[i]=*NNN_Iterator1;
  };
  t.stop();
   
  std::cout << "computed" << std::endl
  << query_point_number << " queries in time " << t.time() <<
  " seconds " << std::endl;

return 0;
};

int main() {
  test_benchmark_nearest_neighbour_L2();
  return 0;
};


