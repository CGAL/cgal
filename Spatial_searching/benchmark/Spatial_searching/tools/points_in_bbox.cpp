#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <iostream>
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;


int main() {

  std::cout.precision(17);

  Point_3 p;
  std::vector<Point_3> points;

  while(std::cin >> p){
    points.push_back(p);
  }

  CGAL::Bbox_3 bb = CGAL::bounding_box(points.begin(), points.end()).bbox();
  double dx = bb.xmax() - bb.xmin();
  double dy = bb.ymax() - bb.ymin();
  double dz = bb.zmax() - bb.zmin();
  
  for(int i = 0; i < points.size(); i++){
    double  rx  = CGAL::get_default_random().get_double();
    double  ry  = CGAL::get_default_random().get_double();
    double  rz  = CGAL::get_default_random().get_double();
    std::cout << bb.xmin() + dx * rx << " ";
    std::cout << bb.ymin() + dy * ry << " ";
    std::cout << bb.zmin() + dz * rz << std::endl;
  } 
  
  return 0;
}

  
