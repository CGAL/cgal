#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian<double> R;
typedef CGAL::Point_2<R> Point;

typedef CGAL::Creator_uniform_2<double,Point> Creator;

typedef CGAL::Plane_separator<double> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Weighted_Minkowski_distance<Point, Point> Distance;
typedef CGAL::Orthogonal_standard_search<Traits, Point, Distance> 
Neighbour_search;

typedef std::vector<Traits::Item> Vector;
typedef std::vector<Point> Query_vector;

int main() {
  
  const int data_point_number=1000;
  
  typedef std::list<Point> point_list;
  point_list data_points;

  // generate random data points  
  CGAL::Random_points_in_square_2<Point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  Traits tr;
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // generate random query points
  const int query_point_number=5;
  CGAL::Random_points_in_square_2<Point,Creator> h( 1.0);
  Query_vector query_points;
  CGAL::copy_n( h, query_point_number, std::back_inserter(query_points));

  Distance tr_dist;

  std::vector<Neighbour_search::Item_with_distance> nearest_neighbour;
  
  for (int i=0; i < query_point_number; i++) { 
     
     Neighbour_search N(d, query_points[i], tr_dist, 1, 0.0);
     N.the_k_neighbours(std::back_inserter(nearest_neighbour));
  }
  
  for (int j=0; j < query_point_number; j++) { 
       std::cout << "q= " << query_points[j] << " ";
       std::cout << "nn= "    << *(nearest_neighbour[j].first) << " ";
       std::cout << " d(q, nn)= "
<< tr_dist.inverse_of_transformed_distance(nearest_neighbour[j].second)  
                 << std::endl;
  } 

  return 0;
};



