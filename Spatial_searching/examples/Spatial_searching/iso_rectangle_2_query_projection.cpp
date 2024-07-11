#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box.h>

#include <CGAL/Search_traits_2.h>

// Point_3 to Point_2 projection on the fly
template<class K>
struct Projection_xy_property_map
{
  typedef typename K::Point_3 key_type;
  typedef typename K::Point_2 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  friend value_type get(Projection_xy_property_map<K>, const key_type& k)
  {
    return value_type(k.x(), k.y());
  }
};

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;

typedef CGAL::Random_points_in_cube_3<Point_3> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;

typedef CGAL::Search_traits_2<K>Traits_base;
typedef CGAL::Search_traits_adapter<Point_3, Projection_xy_property_map<K>, Traits_base> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;

int main()
{
  const int N = 1000;

  std::list<Point_3> points;
  points.push_back(Point_3(0, 0, 0));

  Tree tree;
  Random_points_iterator rpg;
  for(int i = 0; i < N; i++)
    tree.insert(*rpg++);

  std::list<Point_3> result;

  // define 2D range query
  Point_2 p(0.2, 0.2);
  Point_2 q(0.7, 0.7);

  // Searching an exact range
  // using default value 0.0 for epsilon fuzziness parameter
  Fuzzy_iso_box exact_range(p,q);
  tree.search( std::back_inserter( result ), exact_range);
  std::cout << "The points in the box [0.2, 0.7]^2 are: " << std::endl;
  std::copy (result.begin(), result.end(), std::ostream_iterator<Point_3>(std::cout,"\n") );
  std::cout << std::endl;

  result.clear();

  // Searching a fuzzy range
  // using value 0.1 for fuzziness parameter
  Fuzzy_iso_box approximate_range(p, q, 0.1);
  tree.search(std::back_inserter( result ), approximate_range);
  std::cout << "The points in the fuzzy box [[0.1, 0.3], [0.6, 0.8]]^2 are: " << std::endl;
  std::copy (result.begin(), result.end(), std::ostream_iterator<Point_3>(std::cout,"\n") );
  std::cout << std::endl;
  return 0;
}
