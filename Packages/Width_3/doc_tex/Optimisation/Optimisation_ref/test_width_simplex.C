#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <iostream>
#include <vector>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>

typedef leda_integer NrType;
typedef CGAL::Homogeneous<NrType> RepClass;
typedef RepClass::RT RT;
typedef CGAL::Point_3<RepClass> Point;
typedef CGAL::Width_default_traits_3<RepClass> Widthtraits;
typedef CGAL::Width_3<Widthtraits> Width;

int main(int argc, char* argv[]) {
  //*** Simplex: Homogeneous< leda_integer > ***
  std::vector<Point> pointlist;
  Point p(2,0,0,1);
  Point q(0,1,0,1);
  Point r(0,0,1,1);
  Point s(0,0,0,1);
  pointlist.push_back(p);
  pointlist.push_back(q);
  pointlist.push_back(r);
  pointlist.push_back(s);

  // Compute width of simplex
  Width simplex(pointlist.begin(),pointlist.end());
  
  // Output of squared width, width-planes and optimal direction
  RT WNum, WDenom;
  simplex.get_squared_width(WNum,WDenom);
  std::cout<<"Squared Width: "<<WNum<<"/"<<WDenom<<std::endl;
  
  std::cout<<"Direction: "<<simplex.get_build_direction()<<std::endl;

  CGAL::Plane_3<RepClass> e1,e2;
  simplex.get_width_planes(e1,e2);
  std::cout<<"Planes:"<<std::endl;
  std::cout<<"E1: "<<e1<<std::endl<<"E2: "<<e2<<std::endl;

  std::cout<<"Number of optimal solutions: "
	   <<simplex.get_number_of_optimal_solutions()<<std::endl;

  return(0);
}





