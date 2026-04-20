#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_vg.h>
#include <CGAL/IO/write_vg.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/property_map.h>

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
  std::vector<Point_3> points, points2;
  std::vector<std::pair<Sphere_3, std::vector<std::size_t>>> regions, regions2;

  CGAL::IO::read_VG(CGAL::data_file_path("points_3/spheres_point_set_3.vg"), points, std::back_inserter(regions));
  CGAL_assertion(regions.size() == 10);
  CGAL_assertion(regions[0].first.squared_radius() > 0.01625);
  CGAL_assertion(regions[0].first.squared_radius() < 0.01626);
  CGAL_assertion(regions[0].second.size() == 510);
  CGAL_assertion(regions[1].first.squared_radius() > 0.0061);
  CGAL_assertion(regions[1].first.squared_radius() < 0.0062);
  CGAL_assertion(regions[1].second.size() == 134);
  CGAL_assertion(regions[2].first.squared_radius() > 0.0231);
  CGAL_assertion(regions[2].first.squared_radius() < 0.0232);
  CGAL_assertion(regions[2].second.size() == 210);
  CGAL_assertion(regions[3].first.squared_radius() > 0.0053);
  CGAL_assertion(regions[3].first.squared_radius() < 0.0054);
  CGAL_assertion(regions[3].second.size() == 169);
  CGAL_assertion(regions[4].first.squared_radius() > 0.0478);
  CGAL_assertion(regions[4].first.squared_radius() < 0.0479);
  CGAL_assertion(regions[4].second.size() == 1503);
  CGAL::IO::write_VG("spheres_tmp.vg", points, regions,
    CGAL::parameters::point_map(CGAL::Identity_property_map<Point_3>()));
  CGAL::IO::read_VG("spheres_tmp.vg", points2, std::back_inserter(regions2));

  CGAL_assertion(points.size() == points2.size());
  CGAL_assertion(regions.size() == regions2.size());
  for (std::size_t i = 0; i < regions.size(); ++i) {
    CGAL_assertion(regions[i].first == regions2[i].first);
    CGAL_assertion(regions[i].second.size() == regions2[i].second.size());
  }

  struct Cylinder {
    Cylinder() {}
    Cylinder(const Line_3 &axis, FT &radius) : axis(axis), radius(radius) {}

    Line_3 axis;
    FT radius;
  };

  auto serializer = [](auto& cyl, unsigned int& type, std::size_t& num_params) {
    std::stringstream ss;
    ss << cyl.axis << " " << cyl.radius;
    type = 1;
    num_params = 7;
    return ss.str();
    };

  auto constructor = [](unsigned int& type, const std::string &params) {
    Cylinder cyl;
    std::stringstream ss(params);
    ss >> cyl.axis >> cyl.radius;
    return cyl;
    };

  std::vector<Point_3> points3, points4;
  std::vector<std::pair<Cylinder, std::vector<std::size_t>>> regions3, regions4;

  CGAL::IO::read_VG(CGAL::data_file_path("points_3/cylinders_point_set_3.vg"), points3, std::back_inserter(regions3), CGAL::parameters::constructor(constructor));
  CGAL_assertion(regions3.size() == 2);
  CGAL_assertion(regions3[0].first.radius > 0.03835);
  CGAL_assertion(regions3[0].first.radius < 0.03836);
  CGAL_assertion(regions3[0].second.size() == 903);
  CGAL_assertion(regions3[1].first.radius > 0.04345);
  CGAL_assertion(regions3[1].first.radius < 0.04346);
  CGAL_assertion(regions3[1].second.size() == 910);
  CGAL_assertion(points3.size() == 1813);

  CGAL::IO::write_VG("cylinders_tmp.vg", points3, regions3,
    CGAL::parameters::point_map(CGAL::Identity_property_map<Point_3>()).serializer(serializer));
  CGAL::IO::read_VG("cylinders_tmp.vg", points4, std::back_inserter(regions4), CGAL::parameters::constructor(constructor));

  CGAL_assertion(points3.size() == points4.size());
  CGAL_assertion(regions3.size() == regions4.size());
  for (std::size_t i = 0; i < regions3.size(); ++i) {
    // exact comparison works as the file is written with the same precision as the input file
    CGAL_assertion(regions3[i].first.axis == regions4[i].first.axis);
    CGAL_assertion(regions3[i].first.radius == regions4[i].first.radius);
    CGAL_assertion(regions3[i].second.size() == regions4[i].second.size());
  }

  return 0;
}
