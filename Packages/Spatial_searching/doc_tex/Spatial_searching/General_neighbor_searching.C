#include <CGAL/Cartesian_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kernel_d/Iso_box_d.h>
#include <CGAL/Manhattan_distance_rectangle_point.h>
#include <CGAL/General_standard_search.h>

#include <vector> 
#include <iostream>

typedef CGAL::Cartesian_d<double> R;
typedef R::Point_d Point;
typedef Point::R::FT NT;

typedef CGAL::Iso_box_d<R> Box;

typedef CGAL::Kd_tree_traits_point<Point> TreeTraits;
typedef CGAL::Manhattan_distance_rectangle_point<Box, Point> Distance;
typedef CGAL::General_standard_search<TreeTraits, Distance, Box> 
Neighbor_search;
  
int main() {

  const int dim=4;
  const int bucket_size=3;
  
  const int data_point_number=500;
  const int neighbor_number=5;
  
  typedef std::list<Point> Point_list;
  Point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        NT v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]=rnd.get_double(-1.0,1.0);
        Point random_point(dim,v,v+dim);
        data_points.push_front(random_point);
  }
  
  TreeTraits tr(bucket_size, 3.0, false);

  typedef CGAL::Kd_tree<TreeTraits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // define query item
  double p[dim];
  double q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  0.1;
        q[i2]=  0.2;
  }
  Point pp(dim,p,p+dim);
  Point qq(dim,q,q+dim);
  Box query_item(pp,qq);

  std::vector<Neighbor_search::Point_with_distance> neighbors1;

  Distance tr_dist(dim);

  Neighbor_search N1(d, query_item, tr_dist, neighbor_number, 10.0, false);
 
  N1.the_k_neighbors(std::back_inserter(neighbors1)); 
  std::cout << 
  "query item= [0.1,0.2]^4 " << std::endl <<  "5 approximate furthest neighbors are: " 
  << std::endl; 
  for (int i=0; i < neighbor_number; ++i) { 
     std::cout << " d(q,fn)= " << 
     tr_dist.inverse_of_transformed_distance(neighbors1[i].second) << 
     " fn= " << *(neighbors1[i].first) << std::endl; 
  }

  std::vector<Neighbor_search::Point_with_distance> neighbors2;

  Neighbor_search N2(d, query_item, tr_dist, neighbor_number, 0.0, false);
 
  N2.the_k_neighbors(std::back_inserter(neighbors2));

  std::cout << 
  "query item= [0.1,0.2]^4 " << std::endl <<  "5 exact furthest neighbors are: " 
  << std::endl; 
  for (int k=0; k < neighbor_number; ++k) { 
     std::cout << " d(q,fn)= " 
     << tr_dist.inverse_of_transformed_distance(neighbors2[k].second) << 
     " fn= " << *(neighbors2[k].first) << std::endl; 
  }
  
  return 0;
}; 
  
 


