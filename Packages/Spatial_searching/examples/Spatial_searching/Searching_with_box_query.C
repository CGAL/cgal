// file: examples/Spatial_searching/Searching_with_box_query.C

#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>

typedef CGAL::Homogeneous_d<CGAL::MP_Float> Kernel;
typedef Kernel::Point_d Point_d;
typedef Kernel::Iso_box_d Iso_box;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Manhattan_distance_iso_box_point<Kernel> L1_distance;
typedef CGAL::K_neighbor_search<Kernel, L1_distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

int main() {
  const int D = 4;
  const int N = 1000;
  const int K = 5;
  
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit(4, 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  // define query
  double pa[D] = {0.5, 0.5, 0.5, 0.5};
  double qa[D] = {0.6, 0.6, 0.6, 0.6};
  Point_d p(D,pa,pa+D,1000.0);
  Point_d q(D,qa,qa+D,1000.0);
  Iso_box query(p,q);

  K_neighbor_search search(tree, query, K);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;

  for(K_neighbor_search::iterator it = search.begin(); it!= search.end(); it++) { 
    std::cout << " d(q,nn)= " << it->second 
	      << " nn= " << it->first << std::endl; 
  }
  return 0;
}
