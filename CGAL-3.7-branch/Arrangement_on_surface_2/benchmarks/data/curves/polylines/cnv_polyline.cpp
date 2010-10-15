#include <fstream>
#include <iostream>
#include <cassert.>

int main(int argc, char * argv[])
{
  if (argc < 1) {
    std::cerr << "Usage: cnv filename" << std::endl;
  }

  std::cout << "FileFormat( \"AcsBenchmark\", 0, 1 )" << std::endl
            << "BenchmarkName( \"" << argv[1] << "\" )" << std::endl
            << "Classification( \"Arrangement\", \"Polylines\", \"BoundedArcs\", \"Dgn\", \"O\", \"Aug-2007\" )" << std::endl
            << std::endl;
  
  std::ifstream is(argv[1]);
#if 1
  int c;
  is >> c;
  while (c-- && !is.eof()) {
#else
  while (!is.eof()) {
#endif
    std::cout << "Polyline_2(";
    int n;
    is >> n;
    assert(n > 1);
    int x, y, last_x = 0, last_y = 0;
    bool first = true;
    while (n-- > 1 && !is.eof()) {
      is >> x >> y;
      if (first) first = false;
      else {
        if ((x == last_x) && (y == last_y)) continue;
      }
      std::cout << "Point_2(" << x << ", " << y << "), ";
      last_x = x;
      last_y = y;
    }
    is >> x >> y;
    if ((x != last_x) || (y != last_y))
      std::cout << "Point_2(" << x << ", " << y << "))" << std::endl;;
  }
  return 0;
}
