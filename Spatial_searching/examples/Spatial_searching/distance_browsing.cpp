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
  bool operator()(const NN_iterator& it) { return ((*it).first)[0]<0;  }
};

// An iterator that only enumerates dD points with positive x-coordinate
typedef CGAL::Filter_iterator<NN_iterator, X_not_positive> NN_positive_x_iterator;

int main() {

  Tree tree;
  tree.insert(Point_d(0,0));
  tree.insert(Point_d(1,1));
  tree.insert(Point_d(0,1));
  tree.insert(Point_d(10,110));
  tree.insert(Point_d(45,0));
  tree.insert(Point_d(0,2340));
  tree.insert(Point_d(0,30));

  Point_d query(0,0);

  NN_incremental_search NN(tree, query);
  NN_positive_x_iterator it(NN.end(), X_not_positive(), NN.begin()), end(NN.end(), X_not_positive());

  std::cout <<  "The first 5 nearest neighbours with positive x-coord are: " << std::endl;
  for (int j=0; (j < 5)&&(it!=end); ++j,++it)
    std::cout <<   (*it).first << "  at squared distance = " << it->second << std::endl;

  return 0;
}
