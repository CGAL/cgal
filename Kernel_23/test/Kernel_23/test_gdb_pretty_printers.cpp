#include <CGAL/config.h>
#include <CGAL/Simple_cartesian.h>

int main() {
  CGAL::Point_3<CGAL::Simple_cartesian<double>> p(1, 2, 3);
  asm("int3"); // Trigger debugger breakpoint
  return 1;
}
