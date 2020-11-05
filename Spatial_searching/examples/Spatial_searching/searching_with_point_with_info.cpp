#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef boost::tuple<Point_3,int>                           Point_and_int;

typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
  CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
  Traits_base>                                              Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;

int main() {
  const unsigned int K = 5;
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
  Random_points_iterator rpit( 1.0);
  std::vector<Point_3> points;
  std::vector<int>     indices;

  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));
  points.push_back(Point_3(*rpit++));

  indices.push_back(0);
  indices.push_back(1);
  indices.push_back(2);
  indices.push_back(3);
  indices.push_back(4);
  indices.push_back(5);
  indices.push_back(6);

  // Insert number_of_data_points in the tree
  Tree tree(
    boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
    boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )
  );
  Point_3 query(0.0, 0.0, 0.0);
  Distance tr_dist;

  // search K nearest neighbours
  K_neighbor_search search(tree, query, K);
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << " d(q, nearest neighbor)=  "
              << tr_dist.inverse_of_transformed_distance(it->second) << " "
              << boost::get<0>(it->first)<< " " << boost::get<1>(it->first) << std::endl;
  }
  return 0;
}
