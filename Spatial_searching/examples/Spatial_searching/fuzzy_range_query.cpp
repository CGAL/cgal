#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits_d<K> Traits;
typedef CGAL::Random_points_in_cube_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;

int main() {
  const int D = 4;
  const int N = 1000;
  // generator for random data points in the square ( (-1000,-1000), (1000,1000) )
  Random_points_iterator rpit(4, 1000.0);

  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(rpit,N));

  // define range query objects
  double  pcoord[D] = { 300, 300, 300, 300 };
  double  qcoord[D] = { 900.0, 900.0, 900.0, 900.0 };
  Point_d p(D, pcoord, pcoord+D);
  Point_d q(D, qcoord, qcoord+D);
  Fuzzy_sphere fs(p, 700.0, 100.0);
  Fuzzy_iso_box fib(p, q, 100.0);

  std::cout << "points approximately in fuzzy range query" << std::endl;
  std::cout << "with center (300.0, 300.0, 300.0, 300.0)" << std::endl;
  std::cout << "and fuzzy radius <200.0,400.0> are:" << std::endl;
  tree.search(std::ostream_iterator<Point_d>(std::cout, "\n"), fs);

  std::cout << "points approximately in fuzzy range query ";
  std::cout << "[<200,4000>,<800,1000>]]^4 are:" << std::endl;

  tree.search(std::ostream_iterator<Point_d>(std::cout, "\n"), fib);
  return 0;
}
