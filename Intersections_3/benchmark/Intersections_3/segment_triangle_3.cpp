#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
using Triangle_3 = Kernel::Triangle_3;

int main()
{
  Point_3 p(0, 0, 0);
  Point_3 q(10, 0, 0);
  Point_3 r(0, 10, 0);

  Point_3 s(1, 1, -1);
  Point_3 t(1, 1, 1);

  Segment_3 seg(s, t);
  Triangle_3 tri(p, q, r);

  bool result = CGAL::do_intersect(seg, tri);

  if (result) {
    std::cout << "Intersection exists." << std::endl;
  } else {
    std::cout << "No intersection." << std::endl;
  }

  return 0;
}
