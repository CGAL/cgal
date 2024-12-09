// $URL$
// $Id$

#include <CGAL/basic.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include "Point.h"  // defines types Point, Construct_coord_iterator
#include "Distance.h"

#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>

typedef CGAL::Search_traits<double, Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Sliding_midpoint<Traits> Sliding_midpoint;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int main() {

  CGAL::Timer t;
  unsigned int NN_number;
  int bucket_size;
  double Eps;
  std::cout << "Enter number of neighbors to be computed \n" ;
  std::cin >> NN_number;
  std::cout << "Enter approximation factor \n" ;
  std::cin >> Eps;

  std::cout << "Enter bucket size \n" ;
  std::cin >> bucket_size;


  char filename_data_points[80], filename_query_points[80];

  std::cout << "Enter input file name containing query points: \n" ;
  std::cin >> filename_query_points;

  std::cout << "Enter input file name containing data points: \n" ;
  std::cin >> filename_data_points;

  std::ifstream in_data_points, in_query_points;
  int data_point_number;
  int query_point_number;
  int N, N_data_points, N_query_points; // dimension of input data

  in_data_points.open(filename_data_points);
  in_data_points >> N_data_points;
  in_data_points >> data_point_number;


  in_query_points.open(filename_query_points);
  in_query_points >> N_query_points;
  in_query_points >> query_point_number;

  assert(N_data_points==N_query_points);
  N=N_data_points;

  std::cout << "nearest neighbor number = " << NN_number << std::endl;
  std::cout << "approximation factor = " << Eps << std::endl;
  std::cout << "dimension = " << N << std::endl;
  std::cout << "query point number = " << query_point_number << std::endl;
  std::cout << "data point number = " << data_point_number << std::endl;

  typedef std::list<Point> point_list;
  point_list query_points, data_points;

  for (int i = 0; i < query_point_number; i++) {
        std::vector<double> p(N);
        for (int j = 0; j < N; j++) {
          in_query_points >> p[j];
        }
        Point Pnt(p[0],p[1]);
        query_points.push_back(Pnt);
  };

 for (int i = 0; i < data_point_number; i++) {
        std::vector<double> p(N);
        for (int j = 0; j < N; j++) {
          in_data_points >> p[j];
        }
        Point Pnt(p[0],p[1]);
        data_points.push_back(Pnt);
  };

  t.reset();t.start();
  // Insert data points in the tree
  Sliding_midpoint sliding_midpoint(bucket_size);
  Tree tree(data_points.begin(), data_points.end(), sliding_midpoint);
  t.stop();

  data_points.clear();

  std::cout << "created binary search tree containing" << std::endl
  << data_point_number << " points in time "
  << t.time() << std::endl;

  t.reset();t.start();
  for(point_list::iterator it = query_points.begin(); it != query_points.end(); ++it) {
     // Initialize the search structure, and search all NN_number neighbors
     Neighbor_search search(tree, *it, NN_number);
  }
  t.stop();

  std::cout << "time per query is "
  << t.time()/query_point_number << std::endl;
  return 0;
}
