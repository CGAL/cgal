#include <CGAL/basic.h>
#include <CGAL/Cartesian_d.h>  
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Splitting_rules.h>
#include <CGAL/Iso_rectangle_d.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <vector>
#include <iostream>
#include <fstream>

typedef CGAL::Cartesian_d<double> R;
typedef R::Point_d Point;

typedef CGAL::Plane_separator<double> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;

typedef CGAL::Iso_rectangle_d<R> box;	

int main() {

  int dim=3;
  int bucket_size=1;
  CGAL::Timer t;
  
  const int data_point_number=1000000;
  
   
  typedef std::list<Point> point_list;
  point_list data_points,res;
  
  // get data points
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) { 
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        data_points.push_front(Random_point);
  }
  
  /*
  ifstream in;
  int data_point_number;
  
  in.open("data.dat");
  in >> data_point_number;
  
  typedef std::list<Point> point_list;
  point_list data_points,res;

  for (int i = 0; i < data_point_number; i++) {
         
        double v[dim]; 
   	in >> v[0];
        in >> v[1];
        in >> v[2];
        Point P(dim,v,v+dim);
        data_points.push_back(P);
  }; 
  */
  
  Traits tr(bucket_size, CGAL::Split_rules::SLIDING_MIDPOINT, 3.0, true);
  typedef CGAL::Kd_tree<Traits> tree;
  
   
  t.reset(); t.start(); 
  tree d(data_points.begin(), data_points.end(), tr);
  t.stop();
  std::cout << "building time=" << t.time() << std::endl;
  
  
  // define range query
  
  double p[dim];
  double q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  0.2;
        q[i2]=  0.7;
  }
  
  Point P(dim,p,p+dim);
  Point Q(dim,q,q+dim);
  box r(P,Q);

  // Searching the box r
  t.reset();t.start();    
  d.search( std::back_inserter( res ), r);
  t.stop();
  std::cout << "searching time=" << t.time() << std::endl;
  
  std::cout << "Number of the points in the box (0.2,0.2,0.2)-(0.7,0.7,0.7) = " <<
  res.size();
  // std::copy (res.begin(),res.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  

  return 0;
};

