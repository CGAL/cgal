// file: examples/Spatial_searching/Searching_with_box_query.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_d;
typedef CGAL::Search_traits_2<Kernel> Traits;
typedef Traits::Iso_box_d Iso_box;

typedef CGAL::Manhattan_distance_iso_box_point<Traits> L1_distance;
typedef CGAL::K_neighbor_search<Traits, L1_distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

int main() {
  const int D = 2;
  const int N = 1;
  const int K = 1;
  
  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());
 
  Point_d p(0, 0);
  Point_d q(1, 1);
  Iso_box query(p,q);

  K_neighbor_search search(tree, query, K);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;

  for(K_neighbor_search::iterator it = search.begin(); it!= search.end(); it++) { 
    std::cout << " d(q,nn)= " << it->second 
	      << " nn= " << it->first << std::endl; 
  }
  return 0;
}
