// example using nearest_neighbour_iterator for Minkowski_norm with general_distance
// benchmark example using 10000 data points and 2000 query points, bucketsize=10
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
#include <CGAL/Nearest_neighbour_general_distance.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>

int test_benchmark_nearest_neighbour_Minkowski_norm() {

typedef CGAL::Cartesian<double> R;
typedef R::Point_3 Point;
typedef CGAL::Creator_uniform_3<double,Point> Creator;
typedef std::vector<Point> Vector;

typedef CGAL::Kernel_traits<Point>::Kernel K;
typedef K::FT NT;

typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Weighted_Minkowski_distance<Traits> Distance_traits;
typedef CGAL::Nearest_neighbour_general_distance<Traits,
		CGAL::Search_nearest_neighbour, Distance_traits>::iterator 
NNN_Iterator;

  CGAL::Timer t;
  int dim=3;
  int point_number=10000;
  int query_point_number=2000;
  int bucket_size=10;
  NT eps=0.0;
  double the_power=1.3;
  std::vector<double> my_weights(dim);
  my_weights[0]=0.0; my_weights[1]=2.0; my_weights[2]=5;
  
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
  
  // Traits tr(bucket_size, CGAL::SLIDING_MIDPOINT, 3.0, true);
  Traits tr(bucket_size, CGAL::MIDPOINT_OF_MAX_SPREAD, 2.0, false);

  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr, true);
  t.stop();

  std::cout << "created binary search tree containing" << std::endl
  << point_number << " random points in the 3-dim unit cube in time " << t.time() <<
  " seconds " << std::endl;
  d.statistics();

  // end of building binary search tree
  
  std::vector<Traits::Item_with_distance> nearest_neighbours(query_point_number);
  Distance_traits tr_dist(the_power,dim,my_weights);

  t.reset(); t.start();
  for (int i=0; i < query_point_number; i++) {
    // one time iterator
    NNN_Iterator NNN_Iterator1(d,query_points[i],tr_dist,0.0);
    nearest_neighbours[i]=*NNN_Iterator1;
  };
  t.stop();
   
  std::cout << "computed" << std::endl
  << query_point_number << " queries in time " << t.time() <<
  " seconds " << std::endl;
  // brute force approach

  // copy data points from vector to list
  Vector the_data_points;
  the_data_points.reserve(point_number);
  std::copy(data_points.begin(),data_points.end(),std::back_inserter(the_data_points));

  std::vector<int> 
	  nearest_neighbours_brute_force_index(query_point_number);

  t.reset(); t.start();
  for (int i=0; i < query_point_number; i++) {
    // one time iterator
    nearest_neighbours_brute_force_index[i]=0;
    NT distance=
		tr_dist.distance(query_points[i],the_data_points[0]);
    for (int j=1; j < point_number; j++) {
        NT new_distance=tr_dist.distance(query_points[i],
			the_data_points[j]);
		if (new_distance<distance) {
			distance=new_distance;
			nearest_neighbours_brute_force_index[i]=j;
		};
	};
  };
  t.stop();

  std::cout << "computed" << std::endl
  << query_point_number << " queries in time " << t.time() <<
  " seconds using brute force" << std::endl;

  // compare the results
  std::cout << "comparing results" << std::endl;
  for (int i=0; i < query_point_number; i++) {
	if (!(*(nearest_neighbours[i].first)==the_data_points[
		nearest_neighbours_brute_force_index[i]])) {
		assert(
		tr_dist.distance(query_points[i],*(nearest_neighbours[i]).first)==
		tr_dist.distance(query_points[i],
		the_data_points[nearest_neighbours_brute_force_index[i]]));
	};
  };
  std::cout << "all results are fine" << std::endl;
  return 0;
};

int main() {
  test_benchmark_nearest_neighbour_Minkowski_norm();
  return 0;
};


