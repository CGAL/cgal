// file: examples/Spatial_searching/Weighted_Minkowski_distance.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Weighted_Minkowski_distance<TreeTraits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

 
int main() {
  const int D = 2;
  const int N = 1;
  const int K = 1;

  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());
 
  Point_d query(0,0);

  double w[2] = { 1.0, 2.0};
  Distance tr_dist(3.14,D,w, w+D);

  K_neighbor_search search(tree, query, K, 0.0, true, tr_dist);

  for (K_neighbor_search::iterator it = search.begin(); it!=search.end(); ++it){
     std::cout << " d(q,nn)= " << (*it).second 
	       << " nn= " << (*it).first << std::endl; 
  }
  return 0;
}
