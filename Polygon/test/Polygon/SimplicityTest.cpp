#define CGAL_POLYGON_DEBUG 1

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <vector>
#include <fstream>
#include <cstdlib>
#include <cassert>

using std::cout;
using std::endl;

bool TestSimplicity(const char* FileName)
// tests the simplicity of the polygon in the file FileName
{
  typedef CGAL::Exact_rational NT;
  typedef CGAL::Cartesian<NT> K;
  typedef CGAL::Point_2<K> Point;

  cout << "-----------------------------------------------" << endl;
  cout << "-      Testing polygon " << FileName << endl;
  cout << "-----------------------------------------------" << endl;
  std::ifstream from(FileName);
  if (!from) {
    std::cerr << "Error: could not open file " << FileName << endl;
    return false;
  }

  bool answer;            // expected answer (0 indicates false, 1 indicates true)
  int n;                  // number of points
  std::vector<Point> polygon;

  CGAL::set_ascii_mode(from);
  from >> answer >> n;
  cout << "  polygon has " << n << " points" << endl;
  for (int i=0; i<n; i++) {
    Point point;
    from >> point;
    cout << "point " << i << " = " << point << endl;
    polygon.push_back(point);
  }

  bool b = CGAL::is_simple_2(polygon.begin(), polygon.end());

  cout << "(Polygon " << FileName << " is simple) == "
       << (b ? "true" : "false")
       << " (expected result: " << (answer ? "true" : "false") << ")" << endl;
  return (answer == b);
}

void TestDegenerateCases()
{
  typedef CGAL::Cartesian<double> K;
  typedef CGAL::Point_2<K> Point;
  std::vector<Point> polygon;

  polygon.push_back(Point(1,1));
  assert(CGAL::is_simple_2(polygon.begin(), polygon.end()));

  polygon.push_back(Point(1,2));
  assert(CGAL::is_simple_2(polygon.begin(), polygon.end()));
}

int main()
{
  TestDegenerateCases();
  
  bool all_correct = true;
  all_correct &= TestSimplicity("data/simple1.dat");
  all_correct &= TestSimplicity("data/simple2.dat");
  all_correct &= TestSimplicity("data/simple3.dat");
  all_correct &= TestSimplicity("data/simple4.dat");
  all_correct &= TestSimplicity("data/simple5.dat");
  all_correct &= TestSimplicity("data/simple6.dat");
  all_correct &= TestSimplicity("data/simple7.dat");
  all_correct &= TestSimplicity("data/simple8.dat");
  all_correct &= TestSimplicity("data/simple9.dat");
  all_correct &= TestSimplicity("data/simple10.dat");
  all_correct &= TestSimplicity("data/simple11.dat");
  all_correct &= TestSimplicity("data/simple12.dat");

  return all_correct ? 0 : 1;
}

