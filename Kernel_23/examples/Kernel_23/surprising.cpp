#include <iostream>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;

int main()
{
  {
    Point_2 p(0, 0.3), q(1, 0.6), r(2, 0.9);
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");   
  }
  {
    Point_2 p(0, 1.0/3.0), q(1, 2.0/3.0), r(2, 1);
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");   
  }  
  {
    Point_2 p(0,0), q(1, 1), r(2, 2);
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");   
  }  
  return 0;
}
