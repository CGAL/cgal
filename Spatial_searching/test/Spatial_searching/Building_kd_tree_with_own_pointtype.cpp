//file: test/Spatial_searching/Building_kd_tree_with_own_pointtype.C

#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <cassert>
#include "Point.h"
#include "Distance.h"
#include "Point_with_info.h"

typedef CGAL::Random_points_in_cube_3<Point, CGAL::Creator_uniform_3<double,Point> >    Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>                                 N_Random_points_iterator;
typedef CGAL::Search_traits<double, Point, const double*, Construct_coord_iterator>     Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance>                            K_neighbor_search;
//typdefs for Point_with_info
typedef Point_with_info_helper<Point>::type                                             Point_with_info;
typedef Point_property_map<Point>                                                       Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,Traits>                       Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Distance>                         Distance_adapter;
typedef CGAL::Orthogonal_k_neighbor_search<Traits_with_info, Distance_adapter>        K_neighbor_search_with_info;

const unsigned int N = 1000;
const unsigned int K = 5;

template <class K_search>
void run(const std::vector<Point>& points)
{
  // Insert number_of_data_points in the tree
  typename K_search::Tree tree(
    boost::make_transform_iterator(points.begin(),Create_point_with_info<typename K_search::Point_d>()),
    boost::make_transform_iterator(points.end(),Create_point_with_info<typename K_search::Point_d>())
  );

  Point query(0.0, 0.0, 0.0);

  // search K nearest neighbors
  K_search search(tree, query, K);

  // do checking
  double dist = 0;
  std::vector<Point> result;

  for(typename K_search::iterator it = search.begin();
      it != search.end();
      it++){
    result.push_back(get_point(it->first));
    if(CGAL::to_double(it->second) > dist) dist = CGAL::to_double(it->second);
  }

  assert(result.size() == K);
  for(std::vector<Point>::const_iterator it = points.begin();
      it != points.end();
      it++){
    if( std::find(result.begin(), result.end(), *it) == result.end()){
      assert(Distance().transformed_distance(query, *it) >= dist);
    }
  }
  std::cout << "done" << std::endl;
}

int main() {
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
  Random_points_iterator rpit( 1.0);

  std::vector<Point> points(N_Random_points_iterator(rpit,0),
                            N_Random_points_iterator(N));

  run<K_neighbor_search>(points);
  run<K_neighbor_search_with_info>(points);
  return 0;
}
