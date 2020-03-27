#include <CGAL/Search_traits.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include "Point.h"  // defines types Point, Construct_coord_iterator
#include "Distance.h"

typedef CGAL::Creator_uniform_3<double,Point> Point_creator;
typedef CGAL::Random_points_in_cube_3<Point, Point_creator> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Dimension_tag<3> D;
typedef CGAL::Search_traits<double, Point, const double*, Construct_coord_iterator, D> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

int main() {
  const int N = 1000;
  const unsigned int K = 5;
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
  Random_points_iterator rpit( 1.0);

  // Insert number_of_data_points in the tree
  Tree tree(N_Random_points_iterator(rpit,0),
            N_Random_points_iterator(N));

  Point query(0.0, 0.0, 0.0);
  Distance tr_dist;

  // search K nearest neighbours
  K_neighbor_search search(tree, query, K);
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << " d(q, nearest neighbor)=  "
              << tr_dist.inverse_of_transformed_distance(it->second) << std::endl;
  }
  // search K furthest neighbour searching, with eps=0, search_nearest=false
  K_neighbor_search search2(tree, query, K, 0.0, false);

  for(K_neighbor_search::iterator it = search2.begin(); it != search2.end(); it++){
    std::cout << " d(q, furthest neighbor)=  "
              << tr_dist.inverse_of_transformed_distance(it->second) << std::endl;
  }
  return 0;
}
