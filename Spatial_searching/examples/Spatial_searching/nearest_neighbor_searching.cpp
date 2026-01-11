#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <list>
#include <cmath>
#include<cstdlib>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_d;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int main(int argc, char* argv[])
{
  //// 1. Allow user to specify N via command line, default to 1
  const unsigned int N = (argc>1) ? std::atoi(argv[1]) : 1;

  // 2. UX Improvement: Tell the user what is happening
  std::cout<<" generating data and searching for the "<<N<<" nearest neighbor(s)"<< std::endl;

  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());

  // Initialize the search structure, and search all N points
  Point_d query(0,0);
  Neighbor_search search(tree, query, N);

   // report the N nearest neighbors and their distance
  // This should sort all N points by increasing distance from origin
  for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    std::cout << it->first << " "<< std::sqrt(it->second) << std::endl;

  return 0;
}
