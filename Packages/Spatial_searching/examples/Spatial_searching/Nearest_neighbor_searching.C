//file: examples/Spatial_searching/Nearest_neighbor_searching.C

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Kd_tree_traits_point_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
 
typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point_2<R> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef std::list<Neighbor_search::Point_with_distance> Neighbors;

int main() {
  const int N = 1000;
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  // Insert number_of_data_points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  Point query(0,0);
  
  // the container for the result
  Neighbors neighbors;
  
  // Initialize the search structure, and search all N points
  Neighbor_search search(tree, query, N);
  // Perform the search
  search.the_k_neighbors(std::back_inserter(neighbors));

  // report the N nearest neighbors and their distance
  // This should sort all N points by increasing distance from origin
  for(Neighbors::iterator it = neighbors.begin(); it != neighbors.end(); ++it){
    std::cout << it->first << " "<< sqrt(it->second) << std::endl;
  }
  return 0;
}
