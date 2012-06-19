#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel> Traits;
typedef Traits::Point_2 Point_2;

int main(){
  std::vector<Point_2> points;
  std::ifstream input("data/CD500");
  
  double x,y;
  while (input >> x >> y){
    points.push_back(Point_2(x,y,3.));
  }
  
  std::vector<Point_2> ch2;
  
  CGAL::convex_hull_2(points.begin(),points.end(),std::back_inserter(ch2),Traits());
  
  std::cout << ch2.size() << std::endl;
}
