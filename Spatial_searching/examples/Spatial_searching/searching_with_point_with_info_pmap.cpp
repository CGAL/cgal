#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef boost::const_associative_property_map<std::map<std::size_t,Point_3> >           My_point_property_map;

typedef CGAL::Random_points_in_cube_3<Point_3>                                          Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                                                   Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t,My_point_property_map,Traits_base>      Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         Tree;
typedef Tree::Splitter                                                  Splitter;
typedef K_neighbor_search::Distance                                     Distance;

int main() {
  const unsigned int K = 5;
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
  Random_points_iterator rpit( 1.0);
  std::map<std::size_t,Point_3> points;

  points[0]=Point_3(*rpit++);
  points[1]=Point_3(*rpit++);
  points[2]=Point_3(*rpit++);
  points[3]=Point_3(*rpit++);
  points[4]=Point_3(*rpit++);
  points[5]=Point_3(*rpit++);
  points[6]=Point_3(*rpit++);

  My_point_property_map ppmap(points);

  // Insert number_of_data_points in the tree
  Tree tree(
    boost::counting_iterator<std::size_t>(0),
    boost::counting_iterator<std::size_t>(points.size()),
    Splitter(),
    Traits(ppmap)
  );
  Point_3 query(0.0, 0.0, 0.0);
  Distance tr_dist(ppmap);

  // search K nearest neighbours
  K_neighbor_search search(tree, query, K,0,true,tr_dist);
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << " d(q, nearest neighbor)=  "
              << tr_dist.inverse_of_transformed_distance(it->second) << " "
              << points[it->first] << " " << it->first << std::endl;
  }
  return 0;
}
