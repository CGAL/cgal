// file: examples/Spatial_searching/Searching_with_box_query.C

#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_d.h>
#include <CGAL/General_standard_search.h>
#include <CGAL/Manhattan_distance_rectangle_point.h>
#include <iostream>

typedef CGAL::Homogeneous_d<CGAL::MP_Float> R;
typedef R::Point_d Point_d;
typedef R::Iso_box_d Iso_box;
typedef Point_d::R::RT NT;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point_d<R> Traits;
typedef CGAL::Manhattan_distance_rectangle_point<Traits,Iso_box> L1_distance;
typedef CGAL::General_standard_search<Traits, L1_distance> NN_standard_search;
typedef NN_standard_search::Tree Tree;
typedef std::list<NN_standard_search::Point_with_distance> Neighbors;

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

  // define query
  double pa[D] = {0.5, 0.5, 0.5, 0.5};
  double qa[D] = {0.6, 0.6, 0.6, 0.6};
  Point_d p(D,pa,pa+D,1000.0);
  Point_d q(D,qa,qa+D,1000.0);
  Iso_box query(p,q);

  Neighbors neighbors;
  
  NN_standard_search NN(tree, query, K);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  NN.the_k_neighbors(std::back_inserter(neighbors));

  for(Neighbors::iterator it = neighbors.begin(); it!= neighbors.end(); it++) { 
    std::cout << " d(q,nn)= " << it->second 
	      << " nn= " << it->first << std::endl; 
  }
  return 0;
}
