//----------------------------------------------------------------------//
// This is just a small sample application for testing the makefile.
//----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <iostream>

typedef CGAL_Cartesian< double >  R;
typedef CGAL_Point_2< R >         Point;

int main()
{
  Point p(0, 0);
  CGAL_set_ascii_mode(std::cout);
  std::cout << "p = " << p << std::endl;
}

