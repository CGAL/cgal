//----------------------------------------------------------------------//
// This is just a small sample application for testing the makefile.
//----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <iostream.h>

typedef CGAL_Cartesian<double> R;
typedef CGAL_Point_2<R> Point;

int main()
{
  Point p(0,0);
  CGAL_set_ascii_mode(cout);
  cout << "p = " << p << endl;
}

