// file: examples/Spatial_searching/Fuzzy_range_query.C

#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_d.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <iostream>

typedef CGAL::MP_Float NT;
typedef CGAL::Homogeneous_d<NT> R;
typedef R::Point_d Point_d;
typedef R::Iso_box_d Iso_box;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point_d<R> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
typedef CGAL::Fuzzy_iso_box<Traits, Iso_box> Fuzzy_iso_box;

int 
main() {
  const int D = 4;
  const int N = 1000;
  
  // generator for random data points in the square ( (-1000,-1000), (1000,1000) ) 
  Random_points_iterator rpit(4, 1000.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  // define spherical range query object
  NT c[D] = { 300, 300, 300, 300 };
  Point_d center(D,c,c+D, 1.0);
  Fuzzy_sphere fs(center, 700.0, 100.0);

  std::cout << "points approximately in fuzzy range query" << std::endl; 
  std::cout << "with center (300.0, 300.0, 300.0, 300.0)" << std::endl;
  std::cout << "and fuzzy radius <200.0,400.0> are:" << std::endl;
  tree.search(std::ostream_iterator<Point_d>(std::cout, "\n"), fs);
  
  // define rectangular range query object
  NT pa[D] = { -100.0, -100.0, -100.0, -100.0 };
  NT qa[D] = { 900.0, 900.0, 900.0, 900.0 };
  Point_d p(D, pa, pa+D, 1.0);
  Point_d q(D, qa, qa+D, 1.0);
  Fuzzy_iso_box fib(p, q, 100.0);

  std::cout << "points approximately in fuzzy range query ";
  std::cout << "[<-200,0>,<800,1000>]]^4 are:" << std::endl;

  tree.search(std::ostream_iterator<Point_d>(std::cout, "\n"), fib);
 
  return 0;
}

  


