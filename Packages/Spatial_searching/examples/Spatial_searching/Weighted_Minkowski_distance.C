// file: examples/Spatial_searching/Weighted_Minkowski_distance.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
 
typedef CGAL::Cartesian_d<double> K;
typedef K::Point_d Point_d;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point_d<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Weighted_Minkowski_distance<Traits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef   std::list<K_neighbor_search::Point_with_distance> Neighbors;
 
int 
main() {
  const int D=4;
  const int N = 100;
  const int nearest_neighbour_number=10;

  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit(4, 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  double w[4] = { 1.0, 2.0, 3.0, 4.0 };
  Distance tr_dist(3.14,D,w, w+D);


  // define query item
  double q[D] = { 0.5, 0.5, 0.5, 0.5 };
  Point_d query_item(D,q,q+D);

  Neighbors neighbours;
  
  K_neighbor_search NN(tree, query_item, nearest_neighbour_number, 0.0, true, tr_dist);
  NN.the_k_neighbors(std::back_inserter(neighbours));

  for (Neighbors::iterator it = neighbours.begin(); it != neighbours.end();++it) { 
     std::cout << " d(q,nn)= " << (*it).second 
	       << " nn= " << (*it).first << std::endl; 
  }
  return 0;
}
