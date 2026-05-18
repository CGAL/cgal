#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Frechet_distance.h>

#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;

bool load_polyline(const std::string &filename, std::vector<Point> &polyline) {
  polyline.clear();
  std::ifstream ifs(filename);

  if (!ifs)
    return false;

  while (ifs.good()) {
    Point p;
    ifs >> p;
    if (ifs.good())
      polyline.push_back(p);
  }

  return true;
}

int main() {

  std::vector<Point> poly1, poly2;

  if (!load_polyline("poly1.txt", poly1) || !load_polyline("poly2.txt", poly2)) {
    std::cout << "input files could not be loaded" << std::endl;
    return -1;
  }


  std::pair<double, double> res = CGAL::bounded_error_Frechet_distance(poly1, poly2, 0.000001);
 std::cout << "Frechet distance: [" << res.first << ", " << res.second << "]" << std::endl;
  return 0;
}
