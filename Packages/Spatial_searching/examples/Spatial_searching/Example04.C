// Approximate spatial searching: Example04.C
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
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/General_standard_search.h>

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

inline double min_distance_to_queryitem(const Point& p,
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

inline double max_distance_to_queryitem(const Point& p,
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
  
typedef CGAL::Creator_uniform_3<double,Point> Creator;

typedef CGAL::Plane_separator<double> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Orthogonal_standard_search<Traits, Point, Point3D_distance> 
NN_orthogonal_search;
typedef CGAL::General_standard_search<Traits, Point, Point3D_distance> 
NN_general_search;

typedef std::vector<Traits::Item> Vector;
typedef std::vector<Point> Query_vector;

int generate_kd_tree(CGAL::Split_rule_enumeration::Split_rule s) {

  int bucket_size=10;
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  
  Traits tr1(bucket_size, s, 3.0, true);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d1(data_points.begin(), data_points.end(), tr1);

  
  
  std::cout << "created kd tree using extended nodes with splitting rule " 
  << s << " containing "  
  << data_point_number << " points. " << std::endl;
  d1.statistics();

  
  Traits tr2(bucket_size, s, 3.0, false);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d2(data_points.begin(), data_points.end(), tr2);

  std::cout << "created kd tree using no extended nodes with splitting rule " 
  << s << " containing "  
  << data_point_number << " points. " << std::endl;
  d2.statistics();
 
  
  // neighbour searching
  const int query_point_number=5;
  Query_vector query_points;
  CGAL::copy_n( g, query_point_number+1, std::back_inserter(query_points));

  
  Point3D_distance tr_dist;

  // nearest neighbour searching using extended nodes
  std::vector<NN_orthogonal_search::Item_with_distance> nearest_neighbours1;
  // nearest_neighbours1.reserve(query_point_number+1);
  
  
  // nearest neighbour searching using no extended nodes
  std::vector<NN_general_search::Item_with_distance> nearest_neighbours2;
  // nearest_neighbours2.reserve(query_point_number+1);
  
  for (int i=1; i < query_point_number+1; ++i) { 
     NN_orthogonal_search NN1(d1, query_points[i], tr_dist, 1, 0.0);
     std::cout << "neighbour searching statistics using extended nodes: " << std::endl;
     NN1.statistics();
     NN1.the_k_neighbours(std::back_inserter(nearest_neighbours1));
     NN_general_search NN2(d2, query_points[i], tr_dist, 1, 0.0, false);
     std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
     NN2.statistics();
     NN2.the_k_neighbours(std::back_inserter(nearest_neighbours2));
  }
  
  std::cout << "results neighbour searching:" << std::endl;

  for (int i=0; i < query_point_number; ++i) { 
     std::cout << " d(q,nearest neighbour)=" << nearest_neighbours1[i].second << 
     " d(q,furthest neighbour)=" << nearest_neighbours2[i].second << std::endl; 
  } 

  return 0;
};

int main() {
  
  generate_kd_tree(CGAL::Split_rule_enumeration::MEDIAN_OF_MAX_SPREAD); 
  generate_kd_tree(CGAL::Split_rule_enumeration::MEDIAN_OF_RECTANGLE); 
  generate_kd_tree(CGAL::Split_rule_enumeration::MIDPOINT_OF_MAX_SPREAD);
  generate_kd_tree(CGAL::Split_rule_enumeration::MIDPOINT_OF_RECTANGLE);
  generate_kd_tree(CGAL::Split_rule_enumeration::FAIR);
  generate_kd_tree(CGAL::Split_rule_enumeration::SLIDING_MIDPOINT); 
  generate_kd_tree(CGAL::Split_rule_enumeration::SLIDING_FAIR);    

  return 0;
};


