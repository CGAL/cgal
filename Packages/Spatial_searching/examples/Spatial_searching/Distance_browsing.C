//file: examples/Spatial_searching/Distance_browsing.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/iterator.h>
#include <iostream>

typedef double NT;
typedef CGAL::Cartesian_d<NT> K;
typedef K::Point_d Point_d;
typedef CGAL::Random_points_in_iso_box_d<Point_d>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Orthogonal_incremental_neighbor_search<K> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree;

// A functor that returns true, iff the x-coordinate of a dD point is not positive
struct X_not_positive {
  bool operator()(const NN_iterator& it)
  {
    return ((*it).first)[0]<=0;
  }
};

// An iterator that only enumerates dD points with positive x-coordinate
typedef CGAL::Filter_iterator<NN_iterator, X_not_positive> NN_positive_x_iterator;

int 
main() {
  const int D = 4;
  const int N = 20;

  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit(D, 1.0);
  
  // Insert N points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  // define query
  double four[D] = { 0.5, 0.5, 0.5, 0.5 };
  Point_d query(D, four, four+D );

  NN_incremental_search NN(tree, query);

  NN_positive_x_iterator it(NN.begin(), NN.end(), X_not_positive(), NN.begin());
  
  std::cout <<  "The first 5 positive nearest neighbours with positive x-coordinate are: " 
	    << std::endl;

  for (int j=0; j < 5; ++j) { 
      std::cout <<   (*it).first << "  " << (*it).second << std::endl;
    ++it;
   }
   
  return 0;
}
