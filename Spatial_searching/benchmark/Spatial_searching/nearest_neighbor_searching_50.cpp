#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Timer.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_d;
typedef CGAL::Search_traits_3<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


int main() {
  const unsigned int N = 50;
  CGAL::Timer timer;
  timer.start();

  std::vector<Point_d> points, queries;
  Point_d p;

  std::ifstream point_stream("points.xyz");
  while(point_stream >> p){
    points.push_back(p);

  }


  std::ifstream query_stream("queries.xyz");
  while(query_stream >> p ){
    queries.push_back(p);

  }
  timer.stop();
  std::cerr << "reading points took " << timer.time() << " sec." << std::endl;


  timer.reset();
  timer.start();
  Tree tree(points.begin(), points.end());
  tree.build();
  timer.stop();
  std::cerr << "tree construction took " << timer.time() << " sec." << std::endl;


  // Initialize the search structure, and search all N points
  
  double d = 0;
  timer.reset();
  timer.start();
  for(int i = 0; i < queries.size(); i++){
    Neighbor_search search(tree, queries[i], 50, 0);

    // report the N nearest neighbors and their distance
    // This should sort all N points by increasing distance from origin
    for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
      //std::cout << it->first << std::endl;
      d += it->first.x();
    }
  }
  timer.stop();
  std::cerr << d << std::endl;
  std::cerr << queries.size() << " queries in " << timer.time() << " sec." << std::endl;

  return 0;
}
