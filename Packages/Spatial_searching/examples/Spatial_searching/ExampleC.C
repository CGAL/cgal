#include <CGAL/basic.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Iso_rectangle_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Random.h>
#include <CGAL/Splitting_rules.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;
typedef CGAL::Vector_d<R> Vector; // for square root onlys
typedef Point::R::FT NT;

typedef CGAL::Iso_rectangle_d<R> Rectangle;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;

// after CGAL/Kernel/function_objectsHd.h

/*
NT squared_distance(const Point& p, const Point& q) 
{ Vector v = p-q; return v.squared_length(); };
*/  

int main() {

  int bucket_size=1;
  const int dim=4;
  
  
  typedef std::list<Point> point_list;
  point_list data_points;
  const int data_point_number=20;
  
  typedef std::vector<Point> point_vector;
  

  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        double v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]= Rnd.get_double(-1000.0,1000.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  Traits tr(bucket_size, CGAL::Split_rules::SLIDING_FAIR, 3, false);

  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  point_vector points_in_rectangular_range_query;
  point_vector points_in_spherical_range_query;

  point_vector points_in_tree;
  
  d.report_all_points(std::back_inserter(points_in_tree));

  // define center point
  double c[dim];
  for (int i1=0; i1<dim; i1++) {
  	c[i1]=  300.0;
  }
  
  Point C(dim,c,c+dim);

  std::cout << "all points are:" << std::endl;
  
  for (int j1=0; j1 < points_in_tree.size(); ++j1) { 
     std::cout << points_in_tree[j1] << "d(C,p)=" << sqrt(CGAL::squared_distance(points_in_tree[j1],C)) << std::endl; 
  }
  
  
  d.search(std::back_inserter(points_in_spherical_range_query),C,700.0,100.0);

  std::cout << "points approximately in spherical range query are:" << std::endl;
  
  for (int j2=0; j2 < points_in_spherical_range_query.size(); ++j2) { 
     std::cout << points_in_spherical_range_query[j2] << std::endl; 
  }
 
 // define range query
  
  double p[dim];
  double q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  -100.0;
        q[i2]=  900.0;
  }
  
  Point P(dim,p,p+dim);
  Point Q(dim,q,q+dim);

  Rectangle query_rectangle(P,Q);

  d.search(std::back_inserter(points_in_rectangular_range_query),query_rectangle,100.0);

  std::cout << "points approximately in rectangular range query [-100,900]^4 are:" << std::endl;

  
  for (int j3=0; j3 < points_in_rectangular_range_query.size(); ++j3) { 
     std::cout << points_in_rectangular_range_query[j3] << std::endl; 
  }
  

  return 0;
};  

  


