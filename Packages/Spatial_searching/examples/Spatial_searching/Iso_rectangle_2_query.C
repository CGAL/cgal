#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Iso_rectangle_2.h>  
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Splitters.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box_d.h>

#include <vector>
#include <iostream>
#include <fstream>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;

typedef CGAL::Creator_uniform_2<double,Point> Creator;

typedef CGAL::Plane_separator<double> Separator;
typedef CGAL::Kd_tree_traits_point<Point> Traits;

typedef CGAL::Iso_rectangle_2<R> Box;

typedef CGAL::Fuzzy_iso_box_d<Point,Box> Fuzzy_box;	

int main() {

  const int dim=2;
  int bucket_size=1;
  
  const int data_point_number=160;
   
  typedef std::list<Point> point_list;
  point_list data_points, res1, res2;
  
  // generate random data points  
  CGAL::Random_points_in_square_2<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));


  Traits tr(bucket_size, 3.0, true);
  typedef CGAL::Kd_tree<Traits> tree;
  
  tree d(data_points.begin(), data_points.end(), tr);
 
  // define range query
  
  double p[dim];
  double q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  0.2;
        q[i2]=  0.7;
  }
  
  Point P(p[0],p[1]);
  Point Q(q[0],q[1]);
  // Box r(P,Q);

  // Searching an exact range
  // using default value 0.0 for epsilon fuzziness paramater
  // Fuzzy_box exact_range(r); replaced by
  Fuzzy_box exact_range(P,Q);
  d.search( std::back_inserter( res1 ), exact_range);


  std::cout << "The points in the box [0.2,0.7]x[0.2,0.7] are: " << std::endl;
  std::copy (res1.begin(),res1.end(),std::ostream_iterator<Point>(std::cout,"\n") );
  std::cout << std::endl;
  
  // Searching a fuzzy range
  // using value 0.1 for fuzziness paramater
  // Fuzzy_box approximate_range(r,0.1); replaced by
  Fuzzy_box approximate_range(P,Q,0.1);
  d.search( std::back_inserter( res2 ), approximate_range);

  std::cout << "The points in the fuzzy box [<0.1-0.3>,<0.6-0.9>]x[<0.1-0.3>,<0.6-0.9>] are: " << std::endl;
  std::copy (res2.begin(),res2.end(),std::ostream_iterator<Point>(std::cout,"\n") );
  std::cout << std::endl;
  return 0;
};

