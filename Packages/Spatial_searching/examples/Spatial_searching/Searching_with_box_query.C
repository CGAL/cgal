#include <CGAL/basic.h>

#include <vector>
#include <numeric>
#include <cassert>
#include <string>

#include <iostream>
#include <fstream> 

#include <CGAL/Kernel_d/Iso_box_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Random.h>
#include <CGAL/Splitters.h>
#include <CGAL/General_standard_search.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Manhattan_distance_rectangle_point.h>

typedef CGAL::Homogeneous_d<CGAL::MP_Float> R;
typedef R::Point_d Point;
typedef Point::R::RT NT;

typedef CGAL::Iso_box_d<R> Iso_box;
typedef CGAL::Plane_separator<NT> Separator;

typedef CGAL::Kd_tree_traits_point<Point> Traits;
typedef CGAL::Manhattan_distance_rectangle_point<Iso_box,Point> L1_distance;
typedef CGAL::General_standard_search<Traits, L1_distance, Iso_box> 
NN_standard_search;
  

int main() {

  int bucket_size=1;
  const int dim=4;
  
  const int data_point_number=100;
  const int nearest_neighbour_number=10;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) {
	    NT v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim,1.0);
        data_points.push_front(Random_point);
  }
  
  Traits tr(bucket_size, NT(3.0), false);
  L1_distance tr_dist(dim);
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);
  
  // define range query
  NT p[dim];
  NT q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0.5;
        q[i]=  0.6;
  }

  Point P(dim,p,p+dim,1000.0);
  Point Q(dim,q,q+dim,1000.0);

  Iso_box query_rectangle(P,Q);

  std::vector<NN_standard_search::Point_with_distance> nearest_neighbours;
  
  NN_standard_search NN(d, query_rectangle, tr_dist, nearest_neighbour_number, 0.0);
  std::cout << "neighbour searching statistics using no extended nodes: " << std::endl;
  NN.statistics(std::cout);
  NN.the_k_neighbors(std::back_inserter(nearest_neighbours));

  for (int j=0; j < nearest_neighbour_number; ++j) { 
     std::cout << " d(q,nn)= " << nearest_neighbours[j].second << 
     " nn= " << *(nearest_neighbours[j].first) << std::endl; 
  }
  
  return 0;
 
}; 
  
 




