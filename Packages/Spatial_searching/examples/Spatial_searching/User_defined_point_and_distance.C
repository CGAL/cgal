//file: examples/Spatial_searching/User_defined_point_and_distance.C

#include <CGAL/basic.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <iostream>

#include "Point.h"  // defines types Point, Construct_coord_iterator, Distance

typedef CGAL::Random_points_in_cube_3<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree_traits_point<double, Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_standard_search<Traits, Distance> NN_orthogonal_search;
typedef NN_orthogonal_search::Tree Tree;
typedef std::list<NN_orthogonal_search::Point_with_distance> Neighbors;


int
main() {
  const int N = 1000;
  
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  // Insert number_of_data_points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  Point query(0.0,0.0,0.0);
  Distance tr_dist;

  Neighbors neighbors;

  // search 5 nearest neighbours
  NN_orthogonal_search NN1(tree, query, 5);
  NN1.the_k_neighbors(std::back_inserter(neighbors));
  for(Neighbors::iterator it = neighbors.begin(); it != neighbors.end(); it++){
    std::cout << " d(q, nearest neighbor)=  " 
	      << tr_dist.inverse_of_transformed_distance(it->second) << std::endl; 
  }

  neighbors.clear();
  // search 5 furthest neighbour searching, with eps=0, search_nearest=false 
  NN_orthogonal_search NN2(tree, query, 5, 0.0, false);
  NN2.the_k_neighbors(std::back_inserter(neighbors));
  
  for(Neighbors::iterator it = neighbors.begin(); it != neighbors.end(); it++){
    std::cout << " d(q, furthest neighbor)=  " 
	      << tr_dist.inverse_of_transformed_distance(it->second) << std::endl; 
  }
  return 0;
}
