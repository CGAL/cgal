// example using nearest_neighbour_iterator for L2 using standard search with Orthogonal Distance
// benchmark example using 10000 data points and 10000 query points
// bucketsize=1
// both generated with Random_points_in_cube_3<Point_3>
// comparing ASPAS to brute force method 

// #include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 


// #include <CGAL/Cartesian.h>
// #include <CGAL/Point_3.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_L2_standard_search_Minkowski_norm.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
// #include <CGAL/squared_distance_3.h>
// #include <time.h>

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
}; //end of class

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

class Point3D_distance
{
public:

inline double distance(const Point& p1, const Point& p2)
{
	double distx= p1.x()-p2.x();
        double disty= p1.y()-p2.y();
        double distz= p1.z()-p2.z();
        return distx*distx+disty*disty+distz*distz;
}

inline double lower_bound_distance_to_box(const Point& p,
					      const CGAL::Kd_tree_rectangle<double>& b) 
{   double distance(0.0);
    double h;
    h=p.x();
    if (h < b.lower(0)) distance += (b.lower(0)-h)*(b.lower(0)-h);
    if (h > b.upper(0)) distance += (h-b.upper(0))*(h-b.upper(0));
    h=p.y();
    if (h < b.lower(1)) distance += (b.lower(1)-h)*(b.lower(1)-h);
    if (h > b.upper(1)) distance += (h-b.upper(1))*(h-b.upper(1));
    h=p.z();
    if (h < b.lower(2)) distance += (b.lower(2)-h)*(b.lower(2)-h);
    if (h > b.upper(2)) distance += (h-b.upper(2))*(h-b.upper(2));
    return distance;
}

inline double upper_bound_distance_to_box(const Point& p,
					      const CGAL::Kd_tree_rectangle<double>& b) 
{   double distance(0.0);
    double h;
    h=p.x();
    if (h >= (b.lower(0)+b.upper(0))/2.0) 
		  distance += (h-b.lower(0))*(h-b.lower(0)); 
	  else
		  distance += (b.upper(0)-h)*(b.upper(0)-h);
    h=p.y();
    if (h >= (b.lower(1)+b.upper(1))/2.0) 
		  distance += (h-b.lower(1))*(h-b.lower(1)); 
	  else
		  distance += (b.upper(1)-h)*(b.upper(1)-h);
    h=p.z();
    if (h >= (b.lower(2)+b.upper(2))/2.0) 
		  distance += (h-b.lower(2))*(h-b.lower(2)); 
	  else
		  distance += (b.upper(2)-h)*(b.upper(2)-h);
    return distance;
}

inline double new_distance(double& dist, double old_off, double new_off,
			int cutting_dimension)  {
		return dist + new_off*new_off - old_off*old_off;
}

inline double transformed_distance(double d) {
		return d*d;
}

inline double inverse_of_transformed_distance(double d) {
return sqrt(d);
}

}; // end of class

typedef CGAL::Kernel_traits<Point>::Kernel K;
typedef K::FT NT;
// typedef double NT;  causes IRIS compiler problem ??

typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Creator_uniform_3<double,Point> Creator;
typedef CGAL::Nearest_neighbour_L2_standard_search_Minkowski_norm <Traits,CGAL::Search_nearest_neighbour,Point3D_distance> 
Nearest_neighbours_type;

// typedef std::vector<Point> Vector;  causes IRIS compiler problem ??
typedef std::vector<Traits::Item> Vector;

NT The_squared_distance(const Point& P, const Point& Q) {
	NT distx= P.x()-Q.x();
        NT disty= P.y()-Q.y();
        NT distz= P.z()-Q.z();
        return distx*distx+disty*disty+distz*distz;
  };

  int test_benchmark_nearest_neighbour_L2() {

  int dim=3;
  CGAL::Timer t;
  // clock_t ticks_start, ticks_stop;
  // int dim=3; 
  int bucket_size=1;
  NT eps=0.0;

  

  const int data_point_number=  100;
  const int query_point_number= 100;

  std::cout << " bucket_size="
  << bucket_size << " eps=" << eps << std::endl;

  
  // generate 10000 data points
  // Prepare a vector for 10000 points.
    
  // Vector data_points;
  // data_points.reserve(point_number);

  typedef std::list<Point> point_list;
  point_list data_points;
  
  // Create 10000 points within a cube.
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));

  // generate 10000 data points
  // Prepare a vector for 10000 points.
    
  
  Vector query_points;
  query_points.reserve(query_point_number);

  // Create 10000 query points within the same cube.
  CGAL::copy_n( g, query_point_number, std::back_inserter(query_points));

  
  
  
  t.reset(); t.start();
  Traits tr(bucket_size, CGAL::Split_rule_enumeration::SLIDING_MIDPOINT, 3.0, true);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);
  t.stop();

  std::cout << "created binary search tree containing" << std::endl
  << query_point_number << " random points in the 3-dim unit cube in time " 
  << t.time() << std::endl;
  d.statistics();

  // end of building binary search tree
  
  std::vector<CGAL::Nearest_neighbour_L2_standard_search_Minkowski_norm <Traits,
	  CGAL::Search_nearest_neighbour, Point3D_distance>::Item_with_distance> nearest_neighbours;
  nearest_neighbours.reserve(query_point_number);

  t.reset(); t.start();
  
  Point3D_distance tr_dist;
 
  for (int i=0; i < query_point_number; i++) { 
    Nearest_neighbours_type NN(d, query_points[i], tr_dist, 1, 0.0);
    NN.the_k_nearest_neighbours(std::back_inserter(nearest_neighbours));
  };
  
  
  
  t.stop();
   
  std::cout << "computed" << std::endl
  << query_point_number << " queries in time " <<  t.time() <<
  " seconds using ASPAS" << std::endl; 

  /*
  // brute force approach

  // copy data points from list to vector
  
  Vector the_data_points;
  the_data_points.reserve(data_point_number);
  std::copy(data_points.begin(),data_points.end(),
            std::back_inserter(the_data_points));

  std::vector<int> 
	  nearest_neighbours_brute_force_index(query_point_number);

  t.start();
  for (int i=0; i < query_point_number; i++) {
    // one time iterator
    nearest_neighbours_brute_force_index[i]=0;
    NT squared_distance=
		The_squared_distance(query_points[i],the_data_points[0]);
    for (int j=1; j < data_point_number; j++) {
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
  << query_point_number << " queries in time " <<  t.time() <<
  " seconds using brute force" << std::endl;

  // compare the results 
  std::cout << "comparing results" << std::endl;
  for (int i=0; i < query_point_number; i++) {  
	if (!((*(nearest_neighbours[i].first))==the_data_points[
		nearest_neighbours_brute_force_index[i]])) {
		assert(
		(The_squared_distance(query_points[i],
				     (*nearest_neighbours[i].first))==
		The_squared_distance(query_points[i],
		the_data_points[nearest_neighbours_brute_force_index[i]]))); 
	};
  };
  std::cout << "all results are fine" << std::endl;  
  */
  return 0;
};

int main() {
  test_benchmark_nearest_neighbour_L2();

  
  /*
  double dummy;
  std::cout << "Enter input to stop: \n" ;
  std::cin >> dummy;
  */
  
  return 0;
};


