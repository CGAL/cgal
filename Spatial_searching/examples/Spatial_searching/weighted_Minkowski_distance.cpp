#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Random_points_in_square_2<Point_d> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Weighted_Minkowski_distance<TreeTraits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;


int main() {
  const int D = 2;
  const int N = 1000;
  const unsigned int K = 5;


  Random_points_iterator rpit( 1.0);

  Tree tree(N_Random_points_iterator(rpit,0),
	    N_Random_points_iterator(N));

  Point_d query(0,0);

  double w[2] = { 1.0, 2.0};
  Distance tr_dist(3.14,D,w, w+D);

  K_neighbor_search search(tree, query, K, 0.0, true, tr_dist);

  for (K_neighbor_search::iterator it = search.begin(); it!=search.end(); ++it){
    std::cout << "Point " << (*it).first << " at distance = " << (*it).second << std::endl;
  }
  return 0;
}
