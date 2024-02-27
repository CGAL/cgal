#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/internal/implicit_shapes_helper.h>
#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/IO/Image_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Image_3.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <array>
#include <deque>
#include <iostream>

namespace IS = CGAL::Isosurfacing;

template <typename K>
void test_cube()
{
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Mesh = CGAL::Surface_mesh<Point>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Value_function_3<Grid>;
  using Gradients = IS::Finite_difference_gradient_3<K>;
  using Domain = IS::Dual_contouring_domain_3<Grid, Values, Gradients>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::array<std::size_t, 4> >;

  using Mesh = CGAL::Surface_mesh<Point>;

auto implicit_function = [](const Point& q) -> FT
{
  return IS::Shapes::box<K>(Point(-1,-1,-1), Point(1,1,1), q);

  // ---
  auto cyl = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(0,0,1), 0.5, q); };
  auto cube = [](const Point& q) { return IS::Shapes::box<K>(Point(-0.5,-0.5,-0.5), Point(0.5,0.5,0.5), q); };
  auto cyl_and_cube = [&](const Point& q) { return IS::Shapes::shape_union<K>(cyl, cube, q); };

  auto sphere = [](const Point& q) { return IS::Shapes::sphere<K>(Point(0,0,0.5), 1, q); };
  return IS::Shapes::shape_difference<K>(cyl_and_cube, sphere, q);

  // ---
  auto box_1 = [](const Point& q) { return IS::Shapes::box<K>(Point(0,0,0), Point(1,1,1), q); };
  auto box_2 = [](const Point& q) { return IS::Shapes::box<K>(Point(0.5,0.5,0.5), Point(1.5,1.5,1.5), q); };
  return IS::Shapes::shape_union<K>(box_1, box_2, q);

  // ---
  return IS::Shapes::box<K>(Point(0,0,0), Point(1,2,3), q);

  // ---
  return IS::Shapes::torus<K>(Point(0,0,0), Vector(0,0,1), 0.25, 1, q);

  // ---
  return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(1,1,1), 1, q);

  // ---
  auto sphere_1 = [](const Point& q) { return IS::Shapes::sphere<K>(Point(0,0,0), 1, q); };
  auto sphere_2 = [](const Point& q) { return IS::Shapes::sphere<K>(Point(0,0,0.5), 1, q); };
  return IS::Shapes::shape_union<K>(sphere_1, sphere_2, q);
};

  const CGAL::Bbox_3 bbox = {-2., -2., -2., 2., 2., 2.};
  const Vector spacing { 0.05, 0.05, 0.05 };
  const FT isovalue = 0;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring (Implicit function) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };
  Values values {implicit_function , grid };
  Gradients gradients { values, 0.01 * spacing[0] };
  Domain domain { grid, values, gradients };

  Point_range points;
  Polygon_range polygons;
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, polygons, CGAL::parameters::do_not_triangulate_faces(true)); // if you change that, change the array

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << polygons.size() << std::endl;

  assert(points.size() && polygons.size());
  assert(!has_duplicate_points(points, polygons));
  assert(!has_duplicate_polygons(points, polygons));
  assert(!has_isolated_vertices(points, polygons));
  assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons));

  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

  CGAL::IO::write_polygon_mesh("DC_implicit_function.off", mesh, CGAL::parameters::stream_precision(17));

  assert(is_manifold(mesh));
  assert(!has_degenerate_faces(mesh));

  std::cout << "Volume: " << CGAL::Polygon_mesh_processing::volume(mesh) << std::endl;
  assert(CGAL::abs(CGAL::Polygon_mesh_processing::volume(mesh) - FT(8)) < 1e-2);
}

template <typename K>
void test_image()
{
  using FT = typename K::FT;
  using Point = typename K::Point_3;

  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Gradients = IS::Finite_difference_gradient_3<K>;
  using Domain = IS::Dual_contouring_domain_3<Grid, Values, Gradients>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  const std::string fname = CGAL::data_file_path("images/skull_2.9.inr");
  const FT isovalue = 2.9;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring (Image) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    assert(false);
    return;
  }

  // convert image to a Cartesian grid
  Grid grid;
  Values values { grid }; // 'values' keeps a reference to the grid
  if(!IS::IO::read_Image_3(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    assert(false);
    return;
  }

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  // fill up values and gradients
  Gradients gradients { values };
  Domain domain { grid, values, gradients };

  // run dual contouring isosurfacing
  Point_range points;
  Polygon_range polygons;
  IS::dual_contouring(domain, isovalue, points, polygons, CGAL::parameters::do_not_triangulate_faces(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << polygons.size() << std::endl;

  assert(points.size() && polygons.size());
  assert(!has_duplicate_points(points, polygons));
  assert(!has_duplicate_polygons(points, polygons));
  assert(!has_isolated_vertices(points, polygons));
}

int main(int, char**)
{
  test_cube<CGAL::Simple_cartesian<double> >();
  test_image<CGAL::Simple_cartesian<double> >();

  test_cube<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test_image<CGAL::Exact_predicates_inexact_constructions_kernel>();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
