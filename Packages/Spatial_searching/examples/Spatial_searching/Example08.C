// Approximate spatial searching: Example08.C
// Example illustrating for each separate splitting rule 
// building a kd-tree using orthogonal priority search
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
#include <CGAL/Orthogonal_priority_search.h>
#include <CGAL/algorithm.h>

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

typedef CGAL::Kd_tree_rectangle<double> Rectangle;
typedef CGAL::Plane_separator<double> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;

typedef CGAL::Orthogonal_priority_search<Traits, Point, Point3D_distance> 
NN_priority_search;

int test_range_searching(CGAL::Split_rule_enumeration::Split_rule s) {

  std::cout << "test started" << std::endl;

  int bucket_size=1;
  const int dim=3;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) { 
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(v[0],v[1],v[2]);
        data_points.push_front(Random_point);
  }
  
  
  Traits tr(bucket_size, s, 3.0, true);

  Point3D_distance tr_dist;

  std::cout << "constructing tree started" << std::endl;
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);
  std::cout << "constructing tree ready" << std::endl;

  double q[dim];
  q[0]=0.5; q[1]=0.5; q[2]=0.5;
  Point query_item(q[0], q[1], q[2]);

  std::vector<NN_priority_search::Item_with_distance> nearest_neighbours; 
  nearest_neighbours.reserve(nearest_neighbour_number);

  NN_priority_search NN(d, query_item, tr_dist, 0.0, true);

  std::vector<NN_priority_search::Item_with_distance>::iterator 
  it = nearest_neighbours.begin();

  CGAL::copy_n(NN.begin(), nearest_neighbour_number, it);
 
  
  NN.statistics();
  

  for (int i=0; i < nearest_neighbour_number; ++i) { 
     std::cout << " d(q,nn)= " << sqrt(nearest_neighbours[i].second)  <<
     // " nn= " << *(nearest_neighbours[i].first) 
     " nn= " << 
     nearest_neighbours[i].first->x()  << " " <<
     nearest_neighbours[i].first->y()  << " " <<
     nearest_neighbours[i].first->z()  << " " 
     << std::endl; 
  }
  std::cout << "test ready" << std::endl;

  return 0;
}; 
  
int main() {
  
  test_range_searching(CGAL::Split_rule_enumeration::MEDIAN_OF_MAX_SPREAD); 
  test_range_searching(CGAL::Split_rule_enumeration::MEDIAN_OF_RECTANGLE); 
  test_range_searching(CGAL::Split_rule_enumeration::MIDPOINT_OF_MAX_SPREAD);
  test_range_searching(CGAL::Split_rule_enumeration::MIDPOINT_OF_RECTANGLE);
  test_range_searching(CGAL::Split_rule_enumeration::FAIR);
  test_range_searching(CGAL::Split_rule_enumeration::SLIDING_MIDPOINT); 
  test_range_searching(CGAL::Split_rule_enumeration::SLIDING_FAIR);    

  return 0;
};


