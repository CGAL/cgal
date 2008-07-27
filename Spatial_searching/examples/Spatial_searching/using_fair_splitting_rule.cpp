#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <cmath>

typedef CGAL::Simple_cartesian<double> R;
typedef R::Point_2 Point_d;
typedef CGAL::Random_points_in_square_2<Point_d> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<R> Traits;
typedef CGAL::Euclidean_distance<Traits> Distance;
typedef CGAL::Fair<Traits> Fair;
typedef CGAL::Orthogonal_k_neighbor_search<Traits,Distance,Fair> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int main() {
  const unsigned int N = 1000;
  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit( 1.0);

  Fair fair(5); // bucket size=5
  // Insert number_of_data_points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N),
	    fair);

  Point_d query(0,0);

  // Initialize the search structure, and search all N points
  Neighbor_search search(tree, query, N);

  // report the N nearest neighbors and their distance
  // This should sort all N points by increasing distance from origin
  for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
    std::cout << it->first << " "<< std::sqrt(it->second) << std::endl;
  }
  return 0;
}
