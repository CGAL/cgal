#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double>                K;
typedef K::Point_2                                    Point_d;
typedef CGAL::Random_points_in_square_2<Point_d>      Random_points_iterator;
typedef CGAL::Search_traits_2<K>                      Traits;
typedef CGAL::Fuzzy_sphere<Traits>                    Fuzzy_circle;
typedef CGAL::Kd_tree<Traits>                         Tree;

int main() {

  const int N = 30;
  Tree tree;
  Random_points_iterator rpg;
  for(int i = 0; i < N; i++){
    tree.insert(*rpg++);
  }

  // fuzziness = 0

  // Note that a fuzziness of 0 does not imply that we can gather exactly all the
  // points within the disk: even with eps=0, the border is a fuzzy zone.
  Point_d center(0., 0.);
  Fuzzy_circle default_range(center, 0.5);

  std::list<Point_d> result;
  tree.search(std::back_inserter( result ), default_range);

  std::cout << "The points in the fuzzy circle centered at (0., 0.) ";
  std::cout << "with fuzzy radius (0.5, 0.5) are: " << std::endl;
  std::copy (result.begin(),result.end(),std::ostream_iterator<Point_d>(std::cout,"\n") );
  std::cout << std::endl;


  // approximate range searching using value 0.4 for fuzziness parameter
  // We do not write into a list but directly in the outpout stream

  std::cout << "The points in the fuzzy circle centered at (0., 0.) ";
  std::cout << "with fuzzy radius (0.1, 0.9) are: " << std::endl;
  Fuzzy_circle approximate_range(center, 0.5, 0.4);
  tree.search(std::ostream_iterator<Point_d>(std::cout,"\n"), approximate_range);

  return 0;
}
