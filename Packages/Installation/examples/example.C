//----------------------------------------------------------------------//
// This is just a small sample application for testing the makefile.
//----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <iostream>

typedef CGAL::Cartesian< double >  R;
typedef CGAL::Point_2< R >         Point;

int main()
{
  Point p(0, 0);
  CGAL::set_ascii_mode(std::cout);
  std::cout << "p = " << p << std::endl;
}

