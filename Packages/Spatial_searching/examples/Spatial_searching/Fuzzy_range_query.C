#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Fuzzy_sphere_d.h>
#include <CGAL/Fuzzy_iso_box_d.h>
#include <CGAL/Kernel_d/Iso_box_d.h>

#include <vector>
#include <iostream>

typedef CGAL::Homogeneous_d<CGAL::MP_Float> R;
typedef R::Point_d Point;

typedef Point::R::RT NT;

typedef CGAL::Iso_box_d<R> Iso_box;
typedef CGAL::Kd_tree_traits_point<Point> Traits;

typedef CGAL::Fuzzy_sphere_d<Point> Sphere;
typedef CGAL::Fuzzy_iso_box_d<Point, Iso_box> Box;

int main() {

  int bucket_size=1;
  const int dim=4;
  
  typedef std::list<Point> Point_list;
  Point_list data_points;
  const int data_point_number=20;
  
  typedef std::vector<Point> Point_vector;

  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  
  for (int i1=0; i1<data_point_number; i1++) {
        NT v[dim];
        for (int i2=0; i2<dim; i2++) v[i2]= Rnd.get_double(-1000.0,1000.0);
        Point random_point(dim,v,v+dim,1.0);
        data_points.push_front(random_point);
  }
  
  Traits tr(bucket_size, NT(3), false);

  typedef CGAL::Kd_tree<Traits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  Point_vector points_in_rectangular_range_query;
  Point_vector points_in_spherical_range_query;

  // define center point
  NT c[dim];
  for (int i1=0; i1<dim; i1++) {
  	c[i1]=  300.0;
  }
  
  Point center(dim,c,c+dim,1.0);
  Sphere s(center,700.0,100.0);
  d.search(std::back_inserter(points_in_spherical_range_query),s);

  std::cout << "points approximately in fuzzy range query" << std::endl; 
  std::cout << "with center (300.0, 300.0, 300.0, 300.0)" << std::endl;
  std::cout << "and fuzzy radius <200.0,400.0> are:" << std::endl;
  
  unsigned int points_in_spherical_range_query_size=
  points_in_spherical_range_query.size();
  for (unsigned int j2=0; j2 < points_in_spherical_range_query_size; ++j2) { 
     std::cout << points_in_spherical_range_query[j2] << std::endl; 
  }
 
 // define range query
  NT p[dim];
  NT q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  -100.0;
        q[i2]=  900.0;
  }
   
  Point pp(dim,p,p+dim,1.0);
  Point qq(dim,q,q+dim,1.0);

  Box query(pp,qq,100.0);

  d.search(std::back_inserter(points_in_rectangular_range_query),query);

  std::cout << "points approximately in fuzzy range query ";
  std::cout << "[<-200,0>,<800,1000>]]^4 are:" << std::endl;

  unsigned int points_in_rectangular_range_query_size=
               points_in_rectangular_range_query.size();
  for (unsigned int j3=0; j3 < points_in_rectangular_range_query_size; ++j3) { 
     std::cout << points_in_rectangular_range_query[j3] << std::endl; 
  }
  
  return 0;
};  

  


