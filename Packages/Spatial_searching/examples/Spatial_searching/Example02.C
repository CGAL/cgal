// Approximate spatial searching: Example02.C
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
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitting_rules.h>


// create own Point type  
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
}; // end of class Point

  typedef CGAL::Plane_separator<double> Separator;
  typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
  typedef CGAL::Creator_uniform_3<double,Point> Creator;

int generate_kd_tree(CGAL::Split_rule_enumeration::Split_rule s) {

  int bucket_size=10;
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  CGAL::Random_points_in_cube_3<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));

  Traits tr(bucket_size, s, 3.0, true);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  std::cout << "created kd tree using splitting rule " << s << " containing "  
  << data_point_number << " points. " << std::endl;
  d.statistics();
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


