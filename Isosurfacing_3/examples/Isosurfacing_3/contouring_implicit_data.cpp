#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

// using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel, CGAL::Isosurfacing::Do_not_cache_positions>;
using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel, CGAL::Isosurfacing::Cache_positions>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients = CGAL::Isosurfacing::Gradient_function_3<Grid>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

// ---
const FT alpha = 5.01;

auto iwp_value = [](const Point& point)
{
  const FT x = alpha * (point.x() + FT(1.0)) * CGAL_PI;
  const FT y = alpha * (point.y() + FT(1.0)) * CGAL_PI;
  const FT z = alpha * (point.z() + FT(1.0)) * CGAL_PI;
  return cos(x)*cos(y) + cos(y)*cos(z) + cos(z)*cos(x) - cos(x)*cos(y)*cos(z);  // isovalue = 0
};
auto iwp_gradient = [](const Point& point)
{
  const FT x = alpha * (point.x() + FT(1.0)) * CGAL_PI;
  const FT y = alpha * (point.y() + FT(1.0)) * CGAL_PI;
  const FT z = alpha * (point.z() + FT(1.0)) * CGAL_PI;

  const FT gx = CGAL_PI * alpha * sin(x) * (cos(y) * (cos(z) - FT(1.0)) - cos(z));
  const FT gy = CGAL_PI * alpha * sin(y) * (cos(x) * (cos(z) - FT(1.0)) - cos(z));
  const FT gz = CGAL_PI * alpha * sin(z) * (cos(x) * (cos(y) - FT(1.0)) - cos(y));
  return Vector(gx, gy, gz);
};

void run_marching_cubes(const Grid& grid,
                        const FT isovalue)
{
  using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;

  CGAL::Real_timer timer;
  timer.start();

  // fill up values
  Values values { iwp_value, grid };
  Domain domain { grid, values };

  // output containers
  Point_range points;
  Polygon_range triangles;

  // run Marching Cubes
  CGAL::Isosurfacing::marching_cubes(domain, isovalue, points, triangles);

  timer.stop();

  std::cout << "Output #vertices (MC): " << points.size() << std::endl;
  std::cout << "Output #triangles (MC): " << triangles.size() << std::endl;
  std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
  CGAL::IO::write_polygon_soup("marching_cubes_implicit.off", points, triangles);
}

void run_dual_contouring(const Grid& grid,
                         const FT isovalue)
{
  using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  CGAL::Real_timer timer;
  timer.start();

  // fill up values and gradients
  Values values { iwp_value, grid };
  Gradients gradients { iwp_gradient, grid };
  Domain domain { grid, values, gradients };

  // output containers
  Point_range points;
  Polygon_range triangles;

  // run Dual Contouring
  CGAL::Isosurfacing::dual_contouring<CGAL::Sequential_tag>(domain, isovalue, points, triangles,
                                                            CGAL::parameters::do_not_triangulate_faces(true));

  timer.stop();

  std::cout << "Output #vertices (DC): " << points.size() << std::endl;
  std::cout << "Output #triangles (DC): " << triangles.size() << std::endl;
  std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_implicit.off", points, triangles);
}

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.;

  const CGAL::Bbox_3 bbox{-1, -1, -1,  1, 1, 1};
  const FT step = 0.02;
  const Vector spacing { step, step, step };
  Grid grid { bbox, spacing };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  run_marching_cubes(grid, isovalue);

  run_dual_contouring(grid, isovalue);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
