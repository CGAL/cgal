// Copyright (c) 2002 Utrecht University
//
// This file is part of an example program for CGAL. This example
// program may be used, distributed and modified without limitation.
//
//file: examples/Spatial_searching/General_neighbor_searching.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/K_neighbor_search.h>

typedef CGAL::Cartesian_d<double> Kernel;
typedef Kernel::Point_d   Point_d;
typedef Kernel::Iso_box_d Iso_box_d;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Manhattan_distance_iso_box_point<Kernel> Distance;
typedef CGAL::K_neighbor_search<Kernel, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int  main() {
  const int D = 4;
  const int N = 10000;
  const int K = 5;

  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit(D, 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  // Define the query
  double p[D] = {0.1, 0.1, 0.1, 0.1};
  double q[D] = {0.2, 0.2, 0.2, 0.2};
  Point_d pp(D,p,p+D);
  Point_d qq(D,q,q+D);
  Iso_box_d query(pp,qq);

  Distance tr_dist;
  Neighbor_search N1(tree, query, K, 10.0, false); // eps=10.0, nearest=false
  
  std::cout << "query = [0.1,0.2]^4 " << std::endl 
	    <<  K << " approximate furthest neighbors are: " << std::endl; 
  for (Neighbor_search::iterator it = N1.begin();it != N1.end();it++) { 
     std::cout << " d(q,fn)= " << tr_dist.inverse_of_transformed_distance(it->second) 
	       << " fn= " << it->first << std::endl; 
  } 
  return 0;
}
