#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Random.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/General_priority_search.h>
#include <CGAL/L1_distance_rectangle_point.h>

#include <vector>
#include <iostream>

typedef CGAL::Homogeneous<CGAL::MP_Float> R;
typedef CGAL::Point_3<R> Point;
typedef Point::R::FT NT;

typedef CGAL::Iso_cuboid_3<R> Rectangle;
typedef CGAL::Plane_separator<NT> Separator;

typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::L1_distance_rectangle_point<Rectangle,Point> L1_distance;
typedef CGAL::General_priority_search<Traits, Rectangle, L1_distance> 
NN_priority_search;
typedef NN_priority_search::iterator NN_iterator;

template <class InputIterator, class Size, class OutputIterator>
OutputIterator get_elements_in_query( InputIterator first, Size& n,
                       OutputIterator result) {
  
  Size number_of_el_to_compute=n;
  Size count=0;
  while( ((*first).second == NT(0)) && (number_of_el_to_compute > 0)) {
    
    number_of_el_to_compute--;
    count++;
    *result = *first;
    first++;
    result++;
  }
  n=count;
  return result;
}
  
int main() {

  int bucket_size=1;
  const int dim=3;
  
  const int data_point_number=100;
  
  typedef std::list<Point> point_list;
  point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        int v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]= Rnd.get_int(-1000,1000);
        Point Random_point(v[0],v[1],v[2],1000);
        data_points.push_front(Random_point);
  }
  


  Traits tr(bucket_size, CGAL::Split_rule_enumeration::SLIDING_FAIR, 3, false);

  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);
  
  // define range query
  int p[dim];
  int q[dim];
  for (int i=0; i<dim; i++) {
  	p[i]=     0;
        q[i]=  1000;
  }
  
  Point P(p[0],p[1],p[2],1000);
  Point Q(q[0],q[1],q[2],1000);

  Rectangle query_rectangle(P,Q);

  L1_distance tr_dist(dim);
  
  std::vector<NN_priority_search::Item_with_distance> elements_in_query; 
  elements_in_query.reserve(data_point_number);

  NN_priority_search NN(d, query_rectangle, tr_dist, NT(0));

  std::vector<NN_priority_search::Item_with_distance>::iterator 
  it = elements_in_query.begin();

  int n=data_point_number;

  get_elements_in_query(NN.begin(), n, it);
   
  std::cout << 
  "The " << n << " items in the range [0.0,1.0] x [0.0,1.0] x [0.0,1.0] x [0.0,1.0] are: " 
  << std::endl;

  for (int j=0; j < n; ++j) { 
     std::cout <<    
     *(elements_in_query[j].first)
     << std::endl; 
  }
  return 0;
};  

  


