#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_priority_search.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Euclidean_distance.h>

#include <vector>  
#include <iostream>

typedef CGAL::Cartesian_d<double> R;
typedef R::Point_d Point;
typedef Point::R::FT NT;

typedef CGAL::Kd_tree_traits_point<Point> Traits;
typedef CGAL::Euclidean_distance<Point> Distance;

typedef CGAL::Orthogonal_priority_search<Traits> 
NN_priority_search;
typedef NN_priority_search::iterator NN_iterator;

template <class InputIterator, class Size, class OutputIterator>
OutputIterator get_n_positive_elements( InputIterator first, Size& n,
		       Size data_point_number,
                       OutputIterator result) {
  
  Size number_of_el_to_compute=n;
  Size count=0;
  while ((data_point_number>0) && (number_of_el_to_compute)) {
    if ((*(*first).first)[0] > 0) {
    	number_of_el_to_compute--;
    	count++;
    	*result = *first;
    	result++;
    }
    first++;
    data_point_number--;
  }
  n=count;
  return result;
}
  
int main() {

  const int dim=4;
  
  const int data_point_number=20;
  
  typedef std::list<Point> Point_list;
  Point_list data_points;
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
	NT v[dim];
	for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end());

  // define query item
  double q[dim];
  for (int i=0; i<dim; i++) {
  	q[i]=0.5; 
  }
  Point query_item(dim,q,q+dim);

  Distance tr_dist(dim);
  
  std::vector<NN_priority_search::Point_with_distance> elements_in_query; 
  elements_in_query.reserve(data_point_number);

  NN_priority_search NN(d, query_item, tr_dist, 0.0);

  std::vector<NN_priority_search::Point_with_distance>::iterator 
  it = elements_in_query.begin();

  int n=10;
  get_n_positive_elements(NN.begin(), n, data_point_number, it);
   
  std::cout << "query point= 4 0.5 0.5 0.5 0.5 "<< std::endl; 
  std::cout << 
  "The first " << n << " positive nearest neighbours are: " 
  << std::endl;

  for (int j=0; j < n; ++j) { 
     std::cout <<    
     *(elements_in_query[j].first)
     << std::endl; 
  }
  
  return 0;
}; 
  



