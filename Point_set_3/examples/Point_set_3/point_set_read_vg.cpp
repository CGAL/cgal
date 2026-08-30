#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/VG.h>

#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <boost/iterator/function_output_iterator.hpp>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;

using Point_set = CGAL::Point_set_3<Point_3>;
using Line_3 = typename Kernel::Line_3;

using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set>;
using Region_type = CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region_for_point_set<Point_set>;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

int main(int, char**) {
  // Loading example point cloud with some sampled cylinders.
  std::ifstream in(CGAL::data_file_path("points_3/cylinders.ply"));

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Point_set point_set;
  in >> point_set;
  in.close();
  std::cout << "* number of input points: " << point_set.size() << std::endl;
  assert(point_set.has_normal_map()); // input should have normals

  // Parameters to detect cylinders in the point set.
  const std::size_t k = 20;
  const FT          max_distance = FT(1) / FT(10);
  const FT          max_angle = FT(25);
  const std::size_t min_region_size = 20;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(point_set, CGAL::parameters::k_neighbors(k));

  Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_cylinder_fit_region(
    point_set,
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size));

  Region_growing region_growing(
    point_set, neighbor_query, region_type);

  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));

  if (regions.empty())
    return EXIT_FAILURE;

  std::cout << regions.size() << "detected cylinders" << std::endl;

  // Serialization function to write the cylinder parameters into the .vg file.
  // The first parameter takes the cylinder primitive from the region growing,
  // the second parameter is the type of the primitive (here 1 for cylinder),
  // and the third parameter is the number of parameters to write (here 7 for a cylinder: 6 for the axis and 1 for the radius).
  auto serialize = [](auto& cyl, unsigned int& type, std::size_t& num_params) {
    std::stringstream ss;
    ss << cyl.axis << " " << cyl.radius;
    type = 1;
    num_params = 7;
    return ss.str();
    };

  CGAL::IO::write_VG("cylinders_point_set_3.vg", point_set, regions,
    CGAL::parameters::point_map(point_set.point_map()).
    normal_map(point_set.normal_map()).
    serializer(serialize));

  // Cylinder type matching the one used in the serialization function.
  struct Cylinder {
    Cylinder() {}
    Cylinder(const Line_3& axis, FT& radius) : axis(axis), radius(radius) {}

    Line_3 axis;
    FT radius;
  };

  // Constructor function taking the type and a string of parameters to reconstruct the cylinder primitive from the .vg file.
  // The first parameter is not used as only cylinders were exported into the file before.
  auto constructor = [](unsigned int&, const std::string& params) {
    Cylinder cyl;
    std::stringstream ss(params);
    ss >> cyl.axis >> cyl.radius;
    return cyl;
    };

  std::vector<Point_3> points2;
  std::vector<std::pair<Cylinder, std::vector<std::size_t>>> regions2;
  CGAL::IO::read_VG("cylinders_point_set_3.vg", points2, std::back_inserter(regions2), CGAL::parameters::constructor(constructor));

  return EXIT_SUCCESS;
}
