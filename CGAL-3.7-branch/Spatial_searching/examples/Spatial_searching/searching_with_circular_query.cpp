#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Traits> Tree;

int main() {

  const int N=1000;

  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit( 1.0);

  // Insert also the N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
			      N_Random_points_iterator(N));

  // define exact circular range query
  Point center(0.2, 0.2);
  Fuzzy_circle exact_range(center, 0.2);

  std::list<Point> result;
  tree.search(std::back_inserter( result ), exact_range);

  std::cout << "\nPoints in cirle with center " << center << " and radius 0.2" << std::endl;

  std::list<Point>::iterator it;
  for (it=result.begin(); (it != result.end()); ++it) {
    std::cout << *it << std::endl;
  }

  result.clear();
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_circle approximate_range(center, 0.2, 0.1);

  tree.search(std::back_inserter( result ), approximate_range);

  std::cout << "\n\nPoints in cirle with center " << center << " and fuzzy radius [0.1,0.3]" << std::endl;

  for (it=result.begin(); (it != result.end()); ++it) {
    std::cout << *it << std::endl;
  }
  std::cout << "\ndone" << std::endl;
  return 0;
}
