//file: test/Spatial_searching/Building_kd_tree_with_own_pointtype.C

#include <CGAL/basic.h>
#include <CGAL/Search_traits.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <cassert>
#include "Point.h"
#include "Distance.h"

typedef CGAL::Random_points_in_cube_3<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits<double, Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

int main() {
  const unsigned int N = 1000;
  const unsigned int K = 5;
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) ) 
  Random_points_iterator rpit( 1.0);
  
  std::vector<Point> points(N_Random_points_iterator(rpit,0),
			    N_Random_points_iterator(N));
  // Insert number_of_data_points in the tree
  Tree tree(points.begin(), points.end());

  Point query(0.0, 0.0, 0.0);

  // search K nearest neighbours
  K_neighbor_search search(tree, query, K);

  // do checking
  double dist = 0;
  std::vector<Point> result;

  for(K_neighbor_search::iterator it = search.begin();
      it != search.end();
      it++){
    result.push_back(it->first);
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      Distance d;
      assert(d.transformed_distance(query, *it) >= dist);
    }
  }
  std::cout << "done" << std::endl;
  return 0;
}
