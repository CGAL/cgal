// file: examples/Spatial_searching/Distance_browsing.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree;

// A functor that returns true, iff the x-coordinate of a dD point is not positive
struct X_not_positive {
  bool operator()(const NN_iterator& it) { return ((*it).first)[0]<=0;  }
};

// An iterator that only enumerates dD points with positive x-coordinate
typedef CGAL::Filter_iterator<NN_iterator, X_not_positive> NN_positive_x_iterator;

int main() {
  const int D = 2;
  const int N = 1;

  std::list<Point_d> points;
  points.push_back(Point_d(0,0));
  points.push_back(Point_d(1,1));
  points.push_back(Point_d(0,1));
  points.push_back(Point_d(10,110));
  points.push_back(Point_d(45,0));
  points.push_back(Point_d(0,2340));
  points.push_back(Point_d(0,30));

  Tree tree(points.begin(), points.end());
  
  Point_d query(0,0);

  NN_incremental_search NN(tree, query);
  NN_positive_x_iterator it(NN.end(), X_not_positive(), NN.begin());
  
  std::cout <<  "The first 5 nearest neighbours with positive x-coord are: " << std::endl;
  for (int j=0; j < 5; ++j,++it) 
    std::cout <<   (*it).first << "; squared distance =   " << (*it).second << std::endl;
  return 0;
}















