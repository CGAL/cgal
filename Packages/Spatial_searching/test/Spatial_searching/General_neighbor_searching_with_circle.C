#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Euclidean_distance_sphere_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Random.h>
#include <vector>
#include <iostream>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point;
typedef R::FT NT;

typedef R::Circle_2 Circle;

typedef CGAL::Search_traits_2<R> TreeTraits;
typedef CGAL::Euclidean_distance_sphere_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Neighbor_search::Splitter Splitter;
  
int main() {

  const int dim=2;
  const int bucket_size=3;
  
  const int data_point_number=50;
  const int neighbor_number=5;
  
  typedef std::list<Point> Point_list;
  Point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        NT v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]=rnd.get_double(-1.0,1.0);
        Point random_point(v[0],v[1]);
        data_points.push_front(random_point);
  }
  
  Splitter split(bucket_size);
 
  Tree d(data_points.begin(), data_points.end(), split);

  // define query item
  double c[dim];
  for (int i2=0; i2<dim; i2++) {
        c[i2]=  0.2;
  }
  Point cc(c[0],c[1]);
  double squared_radius = 0.04;

  Circle query_item(cc,squared_radius);

  std::vector<Neighbor_search::Point_with_transformed_distance> neighbors1;

  Distance tr_dist;

  Neighbor_search N1(d, query_item, neighbor_number, 10.0, true);
 
  N1.the_k_neighbors(std::back_inserter(neighbors1)); 
  std::cout << 
  "query item= sphere centered at (0.2,0.2) with radius 0.2 " << std::endl <<  
  "5 approximate nearest neighbors are: " 
  << std::endl; 
  for (int i=0; i < neighbor_number; ++i) { 
     std::cout << " d(q,fn)= " << 
     tr_dist.inverse_of_transformed_distance(neighbors1[i].second) << 
     " nn= " << neighbors1[i].first<< std::endl; 
  }

  std::vector<Neighbor_search::Point_with_transformed_distance> neighbors2;

  Neighbor_search N2(d, query_item, neighbor_number, 0.0, true);
 
  N2.the_k_neighbors(std::back_inserter(neighbors2));

  std::cout << 
  "5 exact nearest neighbors are: " 
  << std::endl; 
  for (int k=0; k < neighbor_number; ++k) { 
     std::cout << " d(q,fn)= " 
     << tr_dist.inverse_of_transformed_distance(neighbors2[k].second) << 
     " nn= " << neighbors2[k].first << std::endl; 
  }
  
  return 0;
}; 
  
 


