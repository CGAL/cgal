#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_standard_search.h>

#include <vector>
#include <iostream>
 
typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;

typedef CGAL::Creator_uniform_2<double,Point> Creator;
typedef CGAL::Kd_tree_traits_point<Point> TreeTraits;
typedef CGAL::Orthogonal_standard_search<TreeTraits> Neighbor_search;

typedef std::vector<TreeTraits::Point> Vector;

int main() {
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;

  // generate random data points  
  CGAL::Random_points_in_square_2<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  typedef CGAL::Kd_tree<TreeTraits> Tree;
  Tree d(data_points.begin(), data_points.end());

  // generate random query points
  const int query_point_number=5;
  CGAL::Random_points_in_square_2<Point,Creator> h( 1.0);
  Vector query_points;
  CGAL::copy_n( h, query_point_number, std::back_inserter(query_points));

  std::vector<Neighbor_search::Point_with_distance> the_nearest_neighbors;
  
  for (int i=0; i < query_point_number; i++) { 
     Neighbor_search N(d, query_points[i]); 
     N.the_k_neighbors(std::back_inserter(the_nearest_neighbors));
  }
  
  // report query points q, nearest neighbors and their distance
  for (int j=0; j < query_point_number; j++) { 
       std::cout << "q= " << query_points[j] << " ";
       std::cout << "nn= "    << *(the_nearest_neighbors[j].first) << " ";
       std::cout << " d(q, nn)= "
       << sqrt(the_nearest_neighbors[j].second)  
                 << std::endl;
  } 

  return 0;
}



