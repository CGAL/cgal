/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * example2.C - bench mark
 *    Simple example the CGAL KD-tree module.
 *
 * Written by Sariel Har-Peled 
 *            Iddo Hanniel
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cassert>
#include <list>

#include  <CGAL/kdtree_d.h>
#include  <CGAL/Timer.h>
#include  <CGAL/Random.h>

typedef CGAL::Cartesian<double>           K;
typedef K::Point_3                        Point;
typedef CGAL::Kdtree_interface_3d<Point>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>      kd_tree;
typedef kd_tree::Box                      box;
typedef std::list<Point>                  points_list;
  
int main()
{
  CGAL::Kdtree_d<kd_interface>  tree(3);
  CGAL::Timer t;
  const int dim=3;
  
  // const int data_point_number=1000000;
  const int data_point_number=10000;


  
  
  
  
  typedef std::list<Point> point_list;
  point_list data_points,res;
  
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
  
   
  t.reset();t.start();
  tree.build(data_points);
  t.stop();  
  std::cout << "building time =" << t.time() << std::endl;
     
   
  // Searching the box r
  t.reset();t.start();  
  box r(Point(0.2,0.2,0.2), Point(0.7,0.7,0.7) ,3);
  tree.search( std::back_inserter( res ), r );
  t.stop();
  std::cout << "searching time=" << t.time() << std::endl;
  
  
  std::cout << "Number of the points in the box (0.2,0.2,0.2)-(0.7,0.7,0.7) = " <<
  res.size();
  // std::copy (res.begin(),res.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  return 0;
}
