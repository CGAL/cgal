#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point_3.h>
#include <CGAL/Splitters.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/Fuzzy_iso_box.h>

#include <vector>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Iso_cuboid_3 box;

typedef CGAL::Kd_tree_traits_point_3<K> Traits;
typedef CGAL::Fuzzy_iso_box<Traits,box> Fuzzy_box;	
typedef CGAL::Kd_tree<Traits> Tree;
typedef Tree::Splitter Splitter;


int main() {

  const int dim=3;
  int bucket_size=1;
  CGAL::Timer t;
  
  const int data_point_number=10000;
  // const int data_point_number=1000000;
  
   
  typedef std::list<Point> point_list;
  point_list data_points,res1,res2;
  
  // get data points
  
  // add random points of dimension dim to data_points
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<data_point_number; i1++) { 
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(v[0],v[1],v[2]);
        data_points.push_front(Random_point);
  }
  
  /*
  ifstream in;
  int data_point_number;
  
  in.open("data.dat");
  in >> data_point_number;
  
  typedef std::list<Point> point_list;
  point_list data_points, res1, res2;

  for (int i = 0; i < data_point_number; i++) {
         
        double v[dim]; 
   	in >> v[0];
        in >> v[1];
        in >> v[2];
        Point P(dim,v,v+dim);
        data_points.push_back(P);
  }; 
  */
  
Splitter split(bucket_size, 3.0, true);

   
  t.reset(); t.start(); 
  Tree d(data_points.begin(), data_points.end(), split);
  t.stop();
  std::cout << "building time=" << t.time() << std::endl;
  
  
  // define range query
  
  double p[dim];
  double q[dim];
  for (int i2=0; i2<dim; i2++) {
  	p[i2]=  0.2;
        q[i2]=  0.7;
  }
  
  Point P(p[0],p[1],p[2]);
  Point Q(q[0],q[1],q[2]);

  Fuzzy_box exact_range(P,Q);

  // Searching the box r exactly
  t.reset();t.start();    
  d.search( std::back_inserter( res1 ), exact_range);
  t.stop();
  std::cout << "time exact search=" << t.time() << std::endl;
  
  std::cout << "Number of the points in the box (0.2,0.2,0.2)-(0.7,0.7,0.7) = " <<
  res1.size();
  // std::copy (res1.begin(),res1.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  Fuzzy_box approximate_range(P,Q,0.1);

  // Searching the box r approximately
  t.reset();t.start();    
  d.search( std::back_inserter( res2 ), approximate_range);
  t.stop();
  std::cout << "time approximate search=" << t.time() << std::endl;
  
  std::cout << "Number of the points in the box (0.2,0.2,0.2)-(0.7,0.7,0.7) = " <<
  res2.size();
  // std::copy (res.begin(),res.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;

  return 0;
};

