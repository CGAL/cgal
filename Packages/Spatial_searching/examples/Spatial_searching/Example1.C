// example using nearest_neighbour_iterator for L2
// benchmark example using 10000 data points and 10000 query points
// bucketsize=1
// both generated with Random_points_in_cube_3<Point_3>
// comparing ASPAS to brute force method

#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>

#include <iostream>

//#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_L2.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
// #include <CGAL/squared_distance_3.h>

// typedef CGAL::Cartesian<double> R;
// typedef R::Point_3 Point;

// create own Point type (adapted from example3.C from kdtree and Point_3.h)

class Point
{
public:

  class R
  { 
  public:
    typedef double FT;
  };

private:
  double   vec[ 3 ];
  
public: 
  
  Point()
  { 
    for  ( int ind = 0; ind < 3; ind++ )
      vec[ ind ] = 0;
  }

  Point (double& x, double& y, double& z)
  {
    vec[0]=x;
    vec[1]=y;
    vec[2]=z;
  }

  inline
  int dimension() const
  {
    return  3;
  }
 
  inline
  double x() const
  { 
	return vec[ 0 ];
  }

  inline
  double y() const
  { 
 	return vec[ 1 ];
  }
  
  inline
  double z() const
  { 
	return vec[ 2 ];
  }

  inline
  void set_coord(int k, double x)
  {
    vec[ k ] = x;
  }
  
  inline
  double  & operator[](int k)  
  {
    return  vec[ k ];
  }

  inline
  double  operator[](int k) const
  {
    return  vec[ k ];
  }
};

inline
bool
operator!=(const Point& p, const Point& q)
{
  return ( (p[0] != q[0]) || (p[1] != q[1]) || (p[2] != q[2]) ); 
}

inline
bool
operator==(const Point& p, const Point& q)
{
  return ( (p[0] == q[0]) && (p[1] == q[1]) && (p[2] == q[2]) ) ;
}

typedef CGAL::Creator_uniform_3<double,Point> Creator;
typedef std::vector<Point> Vector;

// typedef CGAL::Kernel_traits<Point>::Kernel K;
// typedef K::FT NT;
typedef double NT;

typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Nearest_neighbour_L2<Traits,
	CGAL::Search_nearest_neighbour>::iterator NNN_Iterator;


NT The_squared_distance(const Point& P, const Point& Q) {
	NT distx= P.x()-Q.x();
        NT disty= P.y()-Q.y();
        NT distz= P.z()-Q.z();
        return distx*distx+disty*disty+distz*distz;
  };

  int test_benchmark_nearest_neighbour_L2() {

  CGAL::Timer t;
  int dim=3;
  int point_number=10000;
  int query_point_number=10000; 
  int bucket_size=1;
  NT eps=0.0;

  std::cout << "test parameters: d=" << dim << " point_number=" 
	    << point_number << std::endl;
  std::cout << "query_point_number=" << query_point_number 
				    << " bucket_size="
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
  query_points.reserve(10000);

  // Create 10000 query points within the same cube.
  CGAL::copy_n( g, query_point_number, std::back_inserter(query_points));

  t.reset(); t.start();
  
  Traits tr(bucket_size, CGAL::SLIDING_MIDPOINT, 3.0, true);
  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr, true);
  t.stop();

  std::cout << "created binary search tree containing" << std::endl
  << point_number << " random points in the 3-dim unit cube in time " 
  << t.time() <<
  " seconds " << std::endl;
  d.statistics();

  // end of building binary search tree
  
  std::vector<Traits::Item_with_distance> 
  nearest_neighbours(query_point_number);

  // CGAL::Timer t1;
  // CGAL::Timer t2;
  t.reset(); t.start();
  // t1.reset(); t2.reset();
  for (int i=0; i < query_point_number; i++) {
    // one time iterator
    // t1.start();t2.start();
    NNN_Iterator NNN_Iterator1(d,query_points[i],0.0);
	// t1.stop();
	// t2.start();
    nearest_neighbours[i]=*NNN_Iterator1;
	// t2.stop();
  };
  t.stop();
   
  
  std::cout << "computed" << std::endl
  << query_point_number << " queries in time " << t.time() <<
  " seconds using ASPAS" << std::endl;

  /*
  std::cout << "computed " << t1.time() <<
  " seconds using time t1" << std::endl;

  std::cout << "computed " << t2.time() <<
  " seconds using time t1+t2" << std::endl;
  */

  // brute force approach

  // copy data points from vector to list
  
  Vector the_data_points;
  the_data_points.reserve(point_number);
  std::copy(data_points.begin(),data_points.end(),
            std::back_inserter(the_data_points));

  std::vector<int> 
	  nearest_neighbours_brute_force_index(query_point_number);

  t.reset(); t.start();
  for (int i=0; i < query_point_number; i++) {
    // one time iterator
    nearest_neighbours_brute_force_index[i]=0;
    NT squared_distance=
		The_squared_distance(query_points[i],the_data_points[0]);
    for (int j=1; j < point_number; j++) {
        NT new_squared_distance=The_squared_distance(query_points[i],
			the_data_points[j]);
		if (new_squared_distance<squared_distance) {
			squared_distance=new_squared_distance;
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
		The_squared_distance(query_points[i],
				     *(nearest_neighbours[i]).first)==
		The_squared_distance(query_points[i],
		the_data_points[nearest_neighbours_brute_force_index[i]]));
	};
  };
  std::cout << "all results are fine" << std::endl; 
  return 0;
};

int main() {
  test_benchmark_nearest_neighbour_L2();

  
  double dummy;
  std::cout << "Enter input to stop: \n" ;
  std::cin >> dummy;
  
  
  return 0;
};


