// file: examples/Spatial_searching/Searching_with_circular_query.C

#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <iostream>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Random_points_in_square_2<Point_d> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;  
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;

int main() {
  const int N = 160;
  
  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());
 
  std::list<Point_d> result;
  
  // define range query
  Point_d center(0.2, 0.2);

  // Searching a circle centered at c with radius 0.2

  // exact range searching using default value 0.0 for fuzziness paramater
  Fuzzy_circle exact_range(center, 0.2);
  tree.search(std::back_inserter( result ), exact_range);
 
  std::cout << "The points in the circle centered at (0.2,0.2) with radius 0.2 are: " << std::endl;
  std::copy (result.begin(),result.end(),std::ostream_iterator<Point_d>(std::cout,"\n") );

  
  // approximate range searching using value 0.1 for fuzziness parameter
  Fuzzy_circle approximate_range(center, 0.2, 0.1);
  tree.search(std::back_inserter( result ), approximate_range);

  std::cout << "\nThe points in the fuzzy circle centered at (0.2,0.2) with fuzzy radius (0.1,0.3) are: " 
	    << std::endl;
  std::copy (result.begin(),result.end(),std::ostream_iterator<Point_d>(std::cout,"\n") );
  return 0;
}
