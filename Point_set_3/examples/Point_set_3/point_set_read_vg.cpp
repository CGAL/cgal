#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_vg.h>
#include <CGAL/Point_set_3.h>

#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Plane_3 = Kernel::Plane_3;
using Sphere_3 = Kernel::Sphere_3;
using Line_3 = typename Kernel::Line_3;
using FT = typename Kernel::FT;

int main(int argc, char** argv) {
  std::vector<Point_3> points;
  std::vector<std::pair<Sphere_3, std::vector<std::size_t>>> regions;

  CGAL::IO::read_VG("spheres_point_set_3.vg", points, std::back_inserter(regions));

  struct Cylinder {
    Cylinder() {}
    Cylinder(const Line_3 &axis, FT &radius) : axis(axis), radius(radius) {}

    Line_3 axis;
    FT radius;
  };

  auto constructor = [](unsigned int& type, const std::string &params) {
    Cylinder cyl;
    std::stringstream ss(params);
    ss >> cyl.axis >> cyl.radius;
    return cyl;
    };

  std::vector<Point_3> points2;
  std::vector<std::pair<Cylinder, std::vector<std::size_t>>> regions2;


  CGAL::IO::read_VG("cylinders_point_set_3.vg", points2, std::back_inserter(regions2), CGAL::parameters::constructor(constructor));
  CGAL_assertion(regions2.size() == 2);
  CGAL_assertion(regions2[0].first.radius > 0.03835);
  CGAL_assertion(regions2[0].first.radius < 0.03836);
  CGAL_assertion(regions2[0].second.size() == 903);
  CGAL_assertion(regions2[1].first.radius > 0.04345);
  CGAL_assertion(regions2[1].first.radius < 0.04346);
  CGAL_assertion(regions2[1].second.size() == 910);
  CGAL_assertion(points2.size() == 1813);

  return 0;
}
