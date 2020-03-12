#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <list>
#include <cmath>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int main() {
  const unsigned int N = 1;

  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());

  Point_d query(0,0);

  // Initialize the search structure, and search all N points

  Neighbor_search search(tree, query, N);

   // report the N nearest neighbors and their distance
  // This should sort all N points by increasing distance from origin
  for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
    std::cout << it->first << " "<< std::sqrt(it->second) << std::endl;
  }


  return 0;
}
