#include <CGAL/basic.h>
#include <vector>
#include <cassert>
#include <iostream>

#include <CGAL/Kd_tree.h>
#include <CGAL/Random.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/algorithm.h>
#include <CGAL/Splitters.h>


#ifdef CUSTOM_POINTS

#include <CGAL/Kd_tree_traits_point.h>
#include "../../examples/Spatial_searching/Point.h" 
typedef CGAL::Kd_tree_traits_point<double, Point, const double*, Construct_coord_iterator> Traits;
#else

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Euclidean_distance<Traits> Distance;
#endif


typedef CGAL::Fair<Traits> Splitter; 
typedef CGAL::Orthogonal_incremental_neighbor_search<Traits, Distance, Splitter>  NN_priority_search;
typedef NN_priority_search::Tree Tree;
typedef NN_priority_search::Splitter Splitter;

int main() {

  std::cout << "test started" << std::endl; 

  int bucket_size=1;
  const int dim=3;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) { 
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(v[0],v[1],v[2]);
        data_points.push_front(Random_point);
  }
  
  Splitter split(bucket_size);

  std::cout << "constructing tree started" << std::endl;
  Tree d(data_points.begin(), data_points.end(), split);
  std::cout << "constructing tree ready" << std::endl;

  double q[dim];
  q[0]=0.5; q[1]=0.5; q[2]=0.5;
  Point query_item(q[0], q[1], q[2]);

  std::vector<NN_priority_search::Point_with_transformed_distance> nearest_neighbours; 
  nearest_neighbours.reserve(nearest_neighbour_number);

  NN_priority_search NN(d, query_item, 0.0, true);

  std::vector<NN_priority_search::Point_with_transformed_distance>::iterator 
  it = nearest_neighbours.begin();

  CGAL::copy_n(NN.begin(), nearest_neighbour_number, it);
 
  
  NN.statistics(std::cout);
  

  for (int i=0; i < nearest_neighbour_number; ++i) { 
     std::cout << " d(q,nn)= " << sqrt(nearest_neighbours[i].second)  <<
     // " nn= " << *(nearest_neighbours[i].first) 
     " nn= " << 
     nearest_neighbours[i].first.x()  << " " <<
     nearest_neighbours[i].first.y()  << " " <<
     nearest_neighbours[i].first.z()  << " " 
     << std::endl; 
  }
  std::cout << "test ready" << std::endl;

  return 0;
}; 
  


