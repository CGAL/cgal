#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Splitters.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere_d.h>

#include <vector>
#include <iostream>
#include <fstream>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;

typedef CGAL::Creator_uniform_2<double,Point> Creator;

typedef CGAL::Plane_separator<double> Separator;
typedef CGAL::Kd_tree_traits_point<Point> Traits;

typedef CGAL::Fuzzy_sphere_d<Point> Fuzzy_circle;

int main() {

  const int dim=2;
  int bucket_size=1;
  
  const int data_point_number=160;
  
   
  typedef std::list<Point> point_list;
  point_list data_points,res1,res2;
  
  // generate random data points  
  CGAL::Random_points_in_square_2<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));

  
  Traits tr(bucket_size, 3.0, true);
  typedef CGAL::Kd_tree<Traits> tree;
  
  tree d(data_points.begin(), data_points.end(), tr);
 
  // define range query
  
  double c[dim];
  for (int i2=0; i2<dim; i2++) {
  	c[i2]=  0.2;
  }
  
  Point center(c[0],c[1]);

  // Searching a circle centered at c with radius 0.2
  // d.search_within_a_radius( std::back_inserter( res ), Center, 0.2);

  // exact range searching using default value 0.0 for fuzziness paramater
  Fuzzy_circle exact_range(center, 0.2);
  d.search(std::back_inserter( res1 ), exact_range);
 
  std::cout << "The points in the circle centered at (0.2,0.2) with radius 0.2 are: " << std::endl;
  std::copy (res1.begin(),res1.end(),std::ostream_iterator<Point>(std::cout,"\n") );
  std::cout << std::endl;
  
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_circle approximate_range(center, 0.2, 0.1);
  d.search(std::back_inserter( res2 ), approximate_range);

  std::cout << "The points in the fuzzy circle centered at (0.2,0.2) with fuzzy radius (0.1,0.3) are: " << std::endl;
  std::copy (res2.begin(),res2.end(),std::ostream_iterator<Point>(std::cout,"\n") );
  std::cout << std::endl;
  
  return 0;
};

