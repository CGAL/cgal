#define CGAL_POLYGON_DEBUG 1

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <CGAL/std/vector>
#include <CGAL/std/fstream>
#include <CGAL/std/cstdlib>
#include <CGAL/std/cassert>

bool TestSimplicity(const char* FileName)
// tests the simplicity of the polygon in the file FileName
{
  typedef CGAL::Cartesian<CGAL::Quotient<CGAL::Gmpz> > R;
  typedef CGAL::Point_2<R> Point;

  CGAL_STD::cout << "-----------------------------------------------" << endl;
  CGAL_STD::cout << "-      Testing polygon " << FileName << endl;
  CGAL_STD::cout << "-----------------------------------------------" << endl;
  CGAL_STD::ifstream from(FileName);
  if (!from) {
    cerr << "Error: could not open file " << FileName << endl;
    return false;
  }

  bool answer;            // expected answer (0 indicates false, 1 indicates true)
  int n;                  // number of points
  CGAL_STD::vector<Point> polygon;

  CGAL::set_ascii_mode(from);
  from >> answer >> n;
  CGAL_STD::cout << "  polygon has " << n << " points" << endl;
  for (int i=0; i<n; i++) {
    Point point;
    from >> point;
    CGAL_STD::cout << "point " << i << " = " << point << endl;
    polygon.push_back(point);
  }

  bool b = CGAL::is_simple_2(polygon.begin(), polygon.end());

  CGAL_STD::cout << "(Polygon " << FileName << " is simple) == "
       << (b ? "true" : "false")
       << " (expected result: " << (answer ? "true" : "false") << ")" << endl;
  return (answer == b);
}

void TestDegenerateCases()
{
  typedef CGAL::Cartesian<double> R;
  typedef CGAL::Point_2<R> Point;
  CGAL_STD::vector<Point> polygon;

  polygon.push_back(Point(1,1));
  assert(CGAL::is_simple_2(polygon.begin(), polygon.end()));

  polygon.push_back(Point(1,2));
  assert(CGAL::is_simple_2(polygon.begin(), polygon.end()));
}

int main()
{
  TestDegenerateCases();

  TestSimplicity("data/simple1.dat");
  TestSimplicity("data/simple2.dat");
  TestSimplicity("data/simple3.dat");
  TestSimplicity("data/simple4.dat");
  TestSimplicity("data/simple5.dat");
  TestSimplicity("data/simple6.dat");
  TestSimplicity("data/simple7.dat");
  TestSimplicity("data/simple8.dat");
  TestSimplicity("data/simple9.dat");
  TestSimplicity("data/simple10.dat");

  return 0;
}

