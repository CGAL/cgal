// How to enter a point in Geomview.
//
//  Sylvain Pion, 2000.

#include <CGAL/Cartesian.h>
#include <iostream>

#if defined(__BORLANDC__) || defined(_MSC_VER)
int main() {
  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
  return 0;
}
#else

#include <CGAL/IO/Geomview_stream.h>

typedef CGAL::Cartesian<double> K;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));

  std::cout << "Please enter a point, by right-clicking on the pickplane"
            << std::endl;

  K::Point_3 p;
  gv >> p;

  std::cout << "Here are the coordinates of the selected point : "
            << p << std::endl;

  return 0;
}
#endif
