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
using Vector_3 = Kernel::Vector_3;

using Point_set = CGAL::Point_set_3<Point_3>;

using Sphere_3 = Kernel::Sphere_3;
using Line_3 = typename Kernel::Line_3;

bool detect_and_export_spheres() {
  using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set>;
  using Region_type = CGAL::Shape_detection::Point_set::Least_squares_sphere_fit_region_for_point_set<Point_set>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

  std::ifstream in(CGAL::data_file_path("points_3/spheres.ply"));

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return false;
  }

  Point_set point_set;
  in >> point_set;
  in.close();
  std::cout << "* number of input points: " << point_set.size() << std::endl;


  // Default parameter values for the data file spheres.ply.
  const std::size_t k = 12;
  const FT          max_distance = FT(1) / FT(100);
  const FT          max_angle = FT(10);
  const std::size_t min_region_size = 50;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(
    point_set, CGAL::parameters::k_neighbors(k));

  Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_sphere_fit_region(
    point_set,
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    point_set, neighbor_query, region_type);

  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));

  if (regions.empty())
    return EXIT_FAILURE;

  std::cout << regions.size() << std::endl;

  CGAL::IO::write_VG("spheres_point_set_3.vg", point_set, regions);

  std::vector<std::pair<Sphere_3, std::vector<std::size_t>>> regions2;
  Point_set point_set2(true);
  CGAL::IO::read_VG("spheres_point_set_3.vg", point_set2, std::back_inserter(regions2));
  CGAL_assertion(regions.size() == regions2.size());

  return true;
}

bool detect_and_export_cylinders() {
  // Similar to the previous function, but with cylinders instead of spheres.
  using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set>;
  using Region_type = CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region_for_point_set<Point_set>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

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
  assert(point_set.size() == 1813);
  assert(point_set.has_normal_map()); // input should have normals

  // Default parameter values for the data file cuble.pwn.
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

  // Create an instance of the region growing class.
  Region_growing region_growing(
    point_set, neighbor_query, region_type);

  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));

  if (regions.empty())
    return EXIT_FAILURE;

  std::cout << regions.size() << std::endl;

  struct Cylinder {
    Cylinder() {}
    Cylinder(const Line_3& axis, FT& radius) : axis(axis), radius(radius) {}

    Line_3 axis;
    FT radius;
  };

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

  auto constructor = [](unsigned int&, const std::string& params) {
    Cylinder cyl;
    std::stringstream ss(params);
    ss >> cyl.axis >> cyl.radius;
    return cyl;
    };

  std::vector<Point_3> points2;
  std::vector<std::pair<Cylinder, std::vector<std::size_t>>> regions2;
  CGAL::IO::read_VG("cylinders_point_set_3.vg", points2, std::back_inserter(regions2), CGAL::parameters::constructor(constructor));
  CGAL_assertion(regions.size() == regions2.size());

  return true;
}

int main(int, char**) {
  if (!detect_and_export_spheres())
    return EXIT_FAILURE;

  if (!detect_and_export_cylinders())
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
