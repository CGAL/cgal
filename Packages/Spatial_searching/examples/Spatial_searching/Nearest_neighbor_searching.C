#include <CGAL/Kd_tree.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;

typedef CGAL::Creator_uniform_2<double,Point> Creator;
typedef CGAL::Kd_tree_traits_point<Point> Traits;
typedef CGAL::Orthogonal_standard_search<Traits> Neighbor_search;

typedef std::vector<Traits::Item> Vector;

int main() {
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;

  // generate random data points  
  CGAL::Random_points_in_square_2<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end());

  // generate random query points
  const int query_point_number=5;
  CGAL::Random_points_in_square_2<Point,Creator> h( 1.0);
  Vector query_points;
  CGAL::copy_n( h, query_point_number, std::back_inserter(query_points));

  std::vector<Neighbor_search::Item_with_distance> nearest_neighbor;
  
  for (int i=0; i < query_point_number; i++) { 
     Neighbor_search N(d, query_points[i]); 
     N.the_k_neighbors(std::back_inserter(nearest_neighbor));
  }
  
  // report query points q, nearest neighbors nn and their distance
  for (int j=0; j < query_point_number; j++) { 
       std::cout << "q= " << query_points[j] << " ";
       std::cout << "nn= "    << *(nearest_neighbor[j].first) << " ";
       std::cout << " d(q, nn)= "
       << sqrt(nearest_neighbor[j].second)  
                 << std::endl;
  } 

  return 0;
}



