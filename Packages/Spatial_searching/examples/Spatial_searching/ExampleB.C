#include <CGAL/basic.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/General_standard_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/basic.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::R::FT NT;

typedef CGAL::Kd_tree_rectangle<NT> Rectangle;
typedef CGAL::Plane_separator<NT> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Weighted_Minkowski_distance<Point, Point> Distance;
typedef CGAL::General_standard_search<Traits, Point, Distance> 
Neighbour_search;
  
int main() {

  const int dim=4;
  const int bucket_size=3;
  const CGAL::Split_rules::Split_rule s=
  CGAL::Split_rules::MEDIAN_OF_MAX_SPREAD;
  
  const int data_point_number=3000;
  const int neighbour_number=5;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        NT v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  Traits tr(bucket_size, s, 3.0, false);

  Distance::Weight_vector w(4);
  w[0]=1.0; w[1]=1.0; w[2]=1.0; w[3]=1.0;

  Distance tr_dist(2,dim,w);

  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // define query item
  double q[dim];
  for (int i=0; i<dim; i++) {
     q[i]=0.0;
  }
  Point query_item(dim,q,q+dim);

  std::vector<Neighbour_search::Item_with_distance> neighbours1;
  neighbours1.reserve(neighbour_number);

  Neighbour_search N1(d, query_item, tr_dist, neighbour_number, 10.0, false);
 
  N1.the_k_neighbours(std::back_inserter(neighbours1));

  std::cout << 
  "query point= 4 0.0 0.0 0.0 0.0" << std::endl <<  "5 approximate furthest neighbours are: " 
  << std::endl; 
  for (int i=0; i < neighbour_number; ++i) { 
     std::cout << " d(q,fn)= " << 
     tr_dist.inverse_of_transformed_distance(neighbours1[i].second) << 
     " fn= " << *(neighbours1[i].first) << std::endl; 
  }

  std::vector<Neighbour_search::Item_with_distance> neighbours2;
  neighbours2.reserve(neighbour_number);

  Neighbour_search N2(d, query_item, tr_dist, neighbour_number, 0.0, false);
 
  N2.the_k_neighbours(std::back_inserter(neighbours2));

  std::cout << 
  "query point= 4 0.0 0.0 0.0 0.0" << std::endl <<  "5 exact furthest neighbours are: " 
  << std::endl; 
  for (int i=0; i < neighbour_number; ++i) { 
     std::cout << " d(q,fn)= " 
     << tr_dist.inverse_of_transformed_distance(neighbours2[i].second) << 
     " fn= " << *(neighbours2[i].first) << std::endl; 
  }
  
  return 0;
}; 
  
 


