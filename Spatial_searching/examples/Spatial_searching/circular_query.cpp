#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Random_points_in_square_2<Point_d> Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Traits> Tree;

int main() {

  const int N=1000;
  Tree tree;
  Random_points_iterator rpg;
  for(int i = 0; i < N; i++){
    tree.insert(*rpg++);
  }

  // define exact circular range query  (fuzziness=0)
  Point_d center(0.2, 0.2);
  Fuzzy_circle exact_range(center, 0.2);

  std::list<Point_d> result;
  tree.search(std::back_inserter( result ), exact_range);

  std::cout << "The points in the circle centered at (0.2,0.2) with radius 0.2 are: " << std::endl;
  std::copy (result.begin(),result.end(),std::ostream_iterator<Point_d>(std::cout,"\n") );
  std::cout << std::endl;

  std::cout << "The points in the fuzzy circle centered at (0.2,0.2) with fuzzy radius (0.1,0.3) are: "
	    << std::endl;
  // approximate range searching using value 0.1 for fuzziness parameter
  // We do not write into a list but directly in the outpout stream
  Fuzzy_circle approximate_range(center, 0.2, 0.1);
  tree.search(std::ostream_iterator<Point_d>(std::cout,"\n"), approximate_range);

  return 0;
}
