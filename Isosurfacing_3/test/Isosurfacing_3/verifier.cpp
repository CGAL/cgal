
#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/IO/polygon_soup_io.h>

/*
Verifier for topological correctness
Based on: Topology Verification for Isosurface Extraction, Tiago Etiene
Data sets: [Marching Cubes cases] [Randomly generated grids] from http://liscustodio.github.io/C_MC33/ 
*/

namespace IS = CGAL::Isosurfacing;

template <typename Grid, typename Values>
void read_iso_volume(const std::string& filename, Grid& grid, Values& values) {
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;
  typename Geom_traits::Construct_point_3 point = grid.geom_traits().construct_point_3_object();
  typename Geom_traits::Construct_iso_cuboid_3 iso_cuboid = grid.geom_traits().construct_iso_cuboid_3_object();
  
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
  }

  // Read the three dimensions
  int nx, ny, nz;
  file.read(reinterpret_cast<char*>(&nx), sizeof(int));
  file.read(reinterpret_cast<char*>(&ny), sizeof(int));
  file.read(reinterpret_cast<char*>(&nz), sizeof(int));

  // Read the bounding box
  // (xmin, xmax, ymin, ymax, zmin, zmax)
  std::array<float, 6> bbox;
  for (int i = 0; i < 6; ++i) {
    file.read(reinterpret_cast<char*>(&bbox[i]), sizeof(float));
  }

  Iso_cuboid_3 span = iso_cuboid(point(bbox[0], bbox[2], bbox[4]),
                                 point(bbox[1], bbox[3], bbox[5]));
  grid = Grid { span, CGAL::make_array<std::size_t>(nx, ny, nz) };

  // Calculate total number of data points
  const std::size_t total_points = nx * ny * nz;

  // Read the volume data (32-bit floats)
  std::vector<float> volume_data(total_points);
  for (std::size_t i = 0; i < total_points; ++i) {
    file.read(reinterpret_cast<char*>(&volume_data[i]), sizeof(float));
    std::cout << volume_data[i] << std::endl;
  }

  for(int x=0; x<nx; ++x)
    for(int y=0; y<ny; ++y)
      for(int z=0; z<nz; ++z)
        values(x, y, z) = volume_data[x + y * nx + z * nx * ny];

  file.close();
}

int read_euler(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
  }
  int euler;
  file >> euler;
  return euler;
}

std::array<int, 2> read_betti(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
  }
  int b0, b1;
  file >> b0;
  file >> b1;
  return {b0, b1};
}

template <typename K>
void verify_euler() {
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Mesh = CGAL::Surface_mesh<Point>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  using Mesh = CGAL::Surface_mesh<Point>;

  const std::size_t num_tests = 10000;

  for (std::size_t i = 0; i < num_tests; i++) {

    Grid grid;
    Values values { grid };
    Domain domain { grid, values };

    read_iso_volume("data/MarchingCubes_cases/Grids/" + std::to_string(i) + "-scalar_field.iso", grid, values);

    Point_range points;
    Polygon_range triangles;
    IS::marching_cubes<CGAL::Sequential_tag>(domain, 0, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

    assert(points.size() && triangles.size());
    assert(!has_duplicate_points(points, triangles));
    assert(!has_duplicate_polygons(points, triangles));
    assert(!has_isolated_vertices(points, triangles));

    assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
    Mesh m;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

    const int euler = euler_characteristic(m);
    const int solution = read_euler("data/MarchingCubes_cases/Grid_invariants/" + std::to_string(i) + "-euler_number.txt");

    if (euler != solution)
      std::cout << "error in test " << i << ": euler " << euler << " != " << solution << std::endl;
  }
}

template <typename K>
void verify_betti() {
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Mesh = CGAL::Surface_mesh<Point>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  using Mesh = CGAL::Surface_mesh<Point>;

  const std::size_t num_tests = 10000;

  for (std::size_t i = 0; i < num_tests; i++) {

    Grid grid;
    Values values { grid };
    Domain domain { grid, values };

    read_iso_volume("data/Closed_Surfaces/Grids/" + std::to_string(i) + "-scalar_field.iso", grid, values);

    Point_range points;
    Polygon_range triangles;
    IS::marching_cubes<CGAL::Sequential_tag>(domain, 0, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

    assert(points.size() && triangles.size());
    assert(!has_duplicate_points(points, triangles));
    assert(!has_duplicate_polygons(points, triangles));
    assert(!has_isolated_vertices(points, triangles));

    assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
    Mesh m;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

    const int b0 = betti_0(m);
    const int b1 = betti_1(m);
    const auto solution = read_betti("data/Closed_Surfaces/InvariantsGrid/" + std::to_string(i) + "-invariant_grid.txt");

    if (b0 != solution[0])
      std::cout << "error in test " << i << ": b0 " << b0 << " != " << solution[0] << std::endl;
    if (b1 != solution[1])
      std::cout << "error in test " << i << ": b1 " << b1 << " != " << solution[1] << std::endl;
  }
}

template <typename K>
void compare_to_reference(const std::string& filename) {
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Mesh = CGAL::Surface_mesh<Point>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  using Mesh = CGAL::Surface_mesh<Point>;

  Grid grid;
  Values values { grid };
  Domain domain { grid, values };

  read_iso_volume(filename, grid, values);

  Point_range points;
  Polygon_range triangles;
  IS::marching_cubes<CGAL::Sequential_tag>(domain, 0, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

  CGAL::IO::write_polygon_soup("verify_tmc.off", points, triangles, CGAL::parameters::stream_precision(17));

  Grid grid_high_res { grid.span().min(), grid.span().max(), std::array<std::size_t, 3>{151, 151, 151} };
  IS::Value_function_3<Grid> values_high_res { values, grid_high_res };
  IS::Marching_cubes_domain_3<Grid, IS::Value_function_3<Grid>> domain_high_res { grid_high_res, values_high_res };

  Point_range points_high_res;
  Polygon_range triangles_high_res;
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain_high_res, 0, points_high_res, triangles_high_res, CGAL::parameters::use_topologically_correct_marching_cubes(false));

  CGAL::IO::write_polygon_soup("verify_reference.off", points_high_res, triangles_high_res);

  write_debug_grid(domain, "verify_cell.off");
}

int main(int, char**)
{
  using K = CGAL::Simple_cartesian<double>;
  using FT = typename K::FT;

  // verify_euler<K>();
  // verify_betti<K>();
  compare_to_reference<K>("data/MarchingCubes_cases/Grids/" + std::to_string(100) + "-scalar_field.iso");
  // compare_to_reference<K>("data/Closed_Surfaces/Grids/" + std::to_string(0) + "-scalar_field.iso");
}