#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/property_map.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Random_points_in_disc_2<Point_2> Random_points;

typedef CGAL::Identity_property_map<Point_2> Point_map;
typedef CGAL::Search_traits_2<Kernel> Search_base;
typedef CGAL::Search_traits_adapter<Point_2, Point_map, Search_base> Search_traits;
typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Search_traits> Tree;

int main ()
{
  Random_points rdpts;
  std::vector<Point_2> pts;
  for (std::size_t i = 0; i < 50; ++ i)
    pts.push_back (*(rdpts ++));

  Tree tree (pts.begin(), pts.end());

  Point_2 center(0., 0.);
  Fuzzy_circle circle (center, 0.5);
  std::vector<Point_2> result;
  tree.search(std::back_inserter(result), circle);
  std::cout << "The points in the fuzzy circle centered at (0., 0.) ";
  std::cout << "with fuzzy radius (0.5) are: " << std::endl;
  for (std::size_t i = 0; i < result.size(); ++ i)
    std::cout << " * " << result[i] << std::endl;

  return 0;
}
