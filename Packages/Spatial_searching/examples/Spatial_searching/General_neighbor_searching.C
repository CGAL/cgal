//file: examples/Spatial_searching/General_neighbor_searching.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree_traits_point_d.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/K_neighbor_search.h>


typedef CGAL::Cartesian_d<double> K;
typedef K::Point_d   Point_d;
typedef K::Iso_box_d Iso_box_d;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point_d<K> TreeTraits;
typedef CGAL::Manhattan_distance_iso_box_point<TreeTraits, Iso_box_d> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef std::list<Neighbor_search::Point_with_distance> Neighbors;

int 
main() {
  const int D = 4;
  const int N = 1000;
  const int K = 5;
  
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit(4, 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  std::cout << "define query" << std::endl;
  double p[D] = {0.1, 0.1, 0.1, 0.1};
  double q[D] = {0.2, 0.2, 0.2, 0.2};
  Point_d pp(D,p,p+D);
  Point_d qq(D,q,q+D);
  Iso_box_d query(pp,qq);

  Neighbors neighbors;
  Neighbors::iterator it;
  Distance tr_dist;

  std::cout << "search" << std::endl;
  Neighbor_search N1(tree, query, K, 10.0, false);
  
  std::cout << "report" << std::endl;
  N1.the_k_neighbors(std::back_inserter(neighbors)); 
 
  std::cout << "query = [0.1,0.2]^4 " << std::endl 
	    <<  K << " approximate furthest neighbors are: " << std::endl; 
  for (it = neighbors.begin();it != neighbors.end();it++) { 
     std::cout << " d(q,fn)= " << tr_dist.inverse_of_transformed_distance(it->second) 
	       << " fn= " << it->first << std::endl; 
  }

  neighbors.clear();

  Neighbor_search N2(tree, query, K, 0.0, false);
 
  N2.the_k_neighbors(std::back_inserter(neighbors));

  std::cout << "query = [0.1,0.2]^4 " << std::endl 
	    <<  K << " exact furthest neighbors are: " << std::endl; 
  for (it = neighbors.begin(); it != neighbors.end(); it++) { 
     std::cout << " d(q,fn)= " << tr_dist.inverse_of_transformed_distance(it->second) 
	       << " fn= " << it->first << std::endl; 
  }  
  return 0;
}
