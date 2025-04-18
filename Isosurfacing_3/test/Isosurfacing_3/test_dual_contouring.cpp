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
#include <CGAL/Random.h>

#include <array>
#include <deque>
#include <iostream>
#include <functional>

namespace IS = CGAL::Isosurfacing;

template <typename K>
void test_implicit_shape(std::function<typename K::FT(const typename K::Point_3&)> implicit_function)
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

  const CGAL::Bbox_3 bbox = {-2., -2., -2., 2., 2., 2.};
  const Vector spacing { 0.05, 0.05, 0.05 };
  const FT isovalue = 0;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring (Implicit function) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };
  Values values { implicit_function, grid };
  Gradients gradients { values, 0.01 * spacing[0] };
  Domain domain { grid, values, gradients };

  Point_range points;
  Polygon_range polygons;
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, polygons, CGAL::parameters::do_not_triangulate_faces(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << polygons.size() << std::endl;

  std::cout << "Checking emptiness..." << std::endl;
  assert(points.size() && polygons.size());
  std::cout << "Checking for duplicate points..." << std::endl;
  assert(!has_duplicate_points(points, polygons));
  std::cout << "Checking for duplicate polygons..." << std::endl;
  assert(!has_duplicate_polygons(points, polygons));
  std::cout << "Checking for isolated vertices..." << std::endl;
  assert(!has_isolated_vertices(points, polygons));
  std::cout << "Checking if the soup is a mesh..." << std::endl;
  assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons));

  std::cout << "Create the mesh..." << std::endl;
  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  std::cout << "Triangulate the mesh..." << std::endl;
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

  std::cout << "Write the polygon mesh..." << std::endl;
  CGAL::IO::write_polygon_mesh("DC_implicit_function.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Check manifoldness..." << std::endl;
  assert(is_manifold(mesh));

  std::cout << "Check for degenerate faces..." << std::endl;
  assert(!has_degenerate_faces(mesh));

  const FT computed_volume = CGAL::Polygon_mesh_processing::volume(mesh);
  std::cout << "Computed Volume: " << computed_volume << std::endl;

  // Brute force estimation of the volume
  std::cout << "Estimating volume using Monte Carlo..." << std::endl;
  const std::size_t num_samples = 1000000;
  CGAL::Random rng(0);

  std::size_t inside_count = 0;
  for (std::size_t i = 0; i < num_samples; ++i) {
    Point sample(rng.get_double(bbox.xmin(), bbox.xmax()),
                 rng.get_double(bbox.ymin(), bbox.ymax()),
                 rng.get_double(bbox.zmin(), bbox.zmax()));
    if (implicit_function(sample) <= isovalue) {
      ++inside_count;
    }
  }

  const FT bbox_volume = (bbox.xmax() - bbox.xmin()) * (bbox.ymax() - bbox.ymin()) * (bbox.zmax() - bbox.zmin());
  const FT monte_carlo_volume = bbox_volume * static_cast<FT>(inside_count) / static_cast<FT>(num_samples);

  std::cout << "Monte Carlo Estimated Volume: " << monte_carlo_volume << std::endl;

  // Compare computed volume with Monte Carlo estimation
  const FT tolerance = 1e-2 * bbox_volume;
  assert(CGAL::abs(computed_volume - monte_carlo_volume) < tolerance);
  std::cout << "Volume check passed!" << std::endl;
}

template <typename K>
void test_box()
{
  auto box_function = [](const typename K::Point_3& q) {
    return IS::Shapes::box<K>(typename K::Point_3(-1, -1, -1), typename K::Point_3(1, 1, 1), q);
  };
  test_implicit_shape<K>(box_function);
}

template <typename K>
void test_union_of_boxes()
{
  auto union_of_boxes_function = [](const typename K::Point_3& q) {
    auto box_1 = [](const typename K::Point_3& q) {
      return IS::Shapes::box<K>(typename K::Point_3(0, 0, 0), typename K::Point_3(1, 1, 1), q);
    };
    auto box_2 = [](const typename K::Point_3& q) {
      return IS::Shapes::box<K>(typename K::Point_3(0.5, 0.5, 0.5), typename K::Point_3(1.5, 1.5, 1.5), q);
    };
    return IS::Shapes::shape_union<K>(box_1, box_2, q);
  };
  test_implicit_shape<K>(union_of_boxes_function);
}

template <typename K>
void test_union_of_boxes_minus_sphere()
{
  auto b1_b2_ms1_function = [](const typename K::Point_3& q) {
    auto box_1 = [](const typename K::Point_3& q) {
      return IS::Shapes::box<K>(typename K::Point_3(0, 0, 0), typename K::Point_3(1, 1, 1), q);
    };
    auto box_2 = [](const typename K::Point_3& q) {
      return IS::Shapes::box<K>(typename K::Point_3(0.5, 0.5, 0.5), typename K::Point_3(1.5, 1.5, 1.5), q);
    };
    auto b1_and_b2 = [&](const typename K::Point_3& q) {
      return IS::Shapes::shape_union<K>(box_1, box_2, q);
    };
    auto s1 = [](const typename K::Point_3& q) {
      return IS::Shapes::sphere<K>(typename K::Point_3(0, 0, 0.5), 1, q);
    };
    return IS::Shapes::shape_difference<K>(b1_and_b2, s1, q);
  };
  test_implicit_shape<K>(b1_b2_ms1_function);
}

template <typename K>
void test_union_of_spheres()
{
  auto union_of_spheres_function = [](const typename K::Point_3& q) {
    auto sphere_1 = [](const typename K::Point_3& q) {
      return IS::Shapes::sphere<K>(typename K::Point_3(0, 0, 0), 1, q);
    };
    auto sphere_2 = [](const typename K::Point_3& q) {
      return IS::Shapes::sphere<K>(typename K::Point_3(0, 0, 0.5), 1, q);
    };
    return IS::Shapes::shape_union<K>(sphere_1, sphere_2, q);
  };
  test_implicit_shape<K>(union_of_spheres_function);
}

template <typename K>
void test_torus()
{
  auto torus_function = [](const typename K::Point_3& q) {
    return IS::Shapes::torus<K>(typename K::Point_3(0, 0, 0), typename K::Vector_3(0, 0, 1), 0.25, 1, q);
  };
  test_implicit_shape<K>(torus_function);
}

// template <typename K>
// void test_infinite_cylinder()
// {
//   auto cylinder_function = [](const typename K::Point_3& q) {
//     return IS::Shapes::infinite_cylinder<K>(typename K::Point_3(0, 0, 0), typename K::Vector_3(1, 1, 1), 1, q);
//   };
//   test_implicit_shape<K>(cylinder_function);
// }

// template <typename K>
// void test_cylinder_and_cube_minus_sphere()
// {
//   auto cyl_cube_function = [](const typename K::Point_3& q) {
//     auto cyl = [](const typename K::Point_3& q) {
//       return IS::Shapes::infinite_cylinder<K>(typename K::Point_3(0, 0, 0), typename K::Vector_3(0, 0, 1), 0.5, q);
//     };
//     auto cube = [](const typename K::Point_3& q) {
//       return IS::Shapes::box<K>(typename K::Point_3(-0.5, -0.5, -0.5), typename K::Point_3(0.5, 0.5, 0.5), q);
//     };
//     auto cyl_and_cube = [&](const typename K::Point_3& q) {
//       return IS::Shapes::shape_union<K>(cyl, cube, q);
//     };
//     auto sphere = [](const typename K::Point_3& q) {
//       return IS::Shapes::sphere<K>(typename K::Point_3(0, 0, 0.5), 1, q);
//     };
//     return IS::Shapes::shape_difference<K>(cyl_and_cube, sphere, q);
//   };
//   test_implicit_shape<K>(cyl_cube_function);
// }

template <typename K>
void test_implicit_shapes()
{
  test_box<K>();
  test_union_of_boxes<K>();
  test_union_of_boxes_minus_sphere<K>();
  test_union_of_spheres<K>();
  test_torus<K>();
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
  if(!IS::IO::convert_image_to_grid(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    assert(false);
    return;
  }

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  const FT step = CGAL::approximate_sqrt(grid.spacing().squared_length()) * 0.01;

  // fill up values and gradients
  Gradients gradients { values, step };
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
  test_implicit_shapes<CGAL::Simple_cartesian<double>>();
  test_image<CGAL::Simple_cartesian<double>>();

  test_implicit_shapes<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test_image<CGAL::Exact_predicates_inexact_constructions_kernel>();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
