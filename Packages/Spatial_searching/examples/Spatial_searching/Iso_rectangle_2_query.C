// file: examples/Spatial_searching/Iso_rectangle_2_query.C

#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box.h>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;	

int 
main() {
  const int N = 1000;
   
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  // Insert number_of_data_points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  std::list<Point> result;
  
  // define range query
  Point p(0.2, 0.2);
  Point q(0.7, 0.7);

  // Searching an exact range
  // using default value 0.0 for epsilon fuzziness paramater
  // Fuzzy_box exact_range(r); replaced by
  Fuzzy_iso_box exact_range(p,q);
  std::cout << "tree.search(..)" << std::endl;
  tree.report_all_points(std::ostream_iterator<Point>(std::cout,"\n"));
  tree.search( std::back_inserter( result ), exact_range);

  std::cout << "The points in the box [0.2,0.7]x[0.2,0.7] are: " << std::endl;
  std::copy (result.begin(), result.end(), std::ostream_iterator<Point>(std::cout,"\n") );
  std::cout << std::endl;
  
  // Searching a fuzzy range
  // using value 0.1 for fuzziness paramater
  std::cout << "The points in the fuzzy box [<0.1-0.3>,<0.6-0.9>]x[<0.1-0.3>,<0.6-0.9>] are: " 
	    << std::endl;
  Fuzzy_iso_box approximate_range(p, q, 0.1);
  tree.search(std::ostream_iterator<Point>(std::cout,"\n") , approximate_range);

  return 0;
}

