// Approximate spatial searching: Example03.C
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
#include <CGAL/Cartesian_d.h>

  typedef CGAL::Cartesian_d<double> R;
  typedef CGAL::Point_d<R> Point;

  typedef CGAL::Plane_separator<double> Separator;
  typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
  
int generate_kd_tree(CGAL::Split_rules::Split_rule s) {

  int bucket_size=10;
  const int dim=100;
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) {
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  Traits tr(bucket_size, s, 3.0, true);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  std::cout << "created kd tree using splitting rule " << s << " containing "  
  << data_point_number << " points. " << std::endl;
  d.statistics();

  return 0;
};

int main() {
  
  generate_kd_tree(CGAL::Split_rules::MEDIAN_OF_MAX_SPREAD); 
  generate_kd_tree(CGAL::Split_rules::MEDIAN_OF_RECTANGLE); 
  generate_kd_tree(CGAL::Split_rules::MIDPOINT_OF_MAX_SPREAD);
  generate_kd_tree(CGAL::Split_rules::MIDPOINT_OF_RECTANGLE);
  generate_kd_tree(CGAL::Split_rules::FAIR);
  generate_kd_tree(CGAL::Split_rules::SLIDING_MIDPOINT); 
  generate_kd_tree(CGAL::Split_rules::SLIDING_FAIR);    

  return 0;
};


