#include <CGAL/Point_3.h>
#include <CGAL/config.h>
#include <CGAL/Simple_cartesian.h>

#include <cstdlib>

int main() {
  CGAL::Point_3<CGAL::Simple_cartesian<double>> p(1, 2, 3);
  std::abort();
  return 1;
}
