#include <array>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;

int main() {
  std::vector<Kernel::Point_3> points;
  std::vector<std::array<std::size_t, 3>> polygons;
  auto pm = CGAL::make_cartesian_converter_property_map<Epick::Point_3>(CGAL::make_identity_property_map(Kernel::Point_3()));
  CGAL::IO::write_polygon_soup("xxx.off", points, polygons, CGAL::parameters::stream_precision(17).point_map(pm));
  return 0;
}
