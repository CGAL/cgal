// file          : test/Spatial_searching/Furthest_neighbor_searching.C

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
 
typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int main() {
  const int N = 1000;
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  Point query(0,0);
  
  // Initialize the search structure, and search all N points
  Neighbor_search search(tree, query, N , 0.0, false);
 
   // report the N nearest neighbors and their distance
  // This should sort all N points by decreasing distance from origin
  double previous_distance=10.0;
  for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
    std::cout << it->first << " "<< sqrt(it->second) << std::endl;
    assert(it->second == CGAL::squared_distance(query,it->first));
    assert(it->second <= previous_distance);
    previous_distance=it->second;
  }
  return 0;
}
