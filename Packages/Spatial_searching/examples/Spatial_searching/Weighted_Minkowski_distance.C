// file: examples/Spatial_searching/Weighted_Minkowski_distance.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
 
typedef CGAL::Cartesian_d<double> Kernel;
typedef Kernel::Point_d Point_d;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Weighted_Minkowski_distance<Kernel> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<Kernel, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

 
int main() {
  const int D = 4;
  const int N = 100;
  const int K = 10;

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

  K_neighbor_search search(tree, query_item, K, 0.0, true, tr_dist);

  for (K_neighbor_search::iterator it = search.begin(); it!=search.end(); ++it){
     std::cout << " d(q,nn)= " << (*it).second 
	       << " nn= " << (*it).first << std::endl; 
  }
  return 0;
}
