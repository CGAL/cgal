#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <cmath>
#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients = CGAL::Isosurfacing::Gradient_function_3<Grid>;
using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

using Mesh = CGAL::Surface_mesh<Point>;

// "Devil" - https://www-sop.inria.fr/galaad/surface/
auto devil_value = [](const Point& point)
{
  const FT x = point.x(), y = point.y(), z = point.z();
  return x*x*x*x + 2*x*x*z*z - 0.36*x*x - y*y*y*y + 0.25*y*y + z*z*z*z;
};

auto devil_gradient = [](const Point& point)
{
  const FT x = point.x(), y = point.y(), z = point.z();

  const FT gx = 4*x*x*x + 4*x*z*z - 0.72*x;
  const FT gy = -4*y*y*y + 0.5*y;
  const FT gz = 4*x*x*z + 4*z*z*z;
  Vector g(gx, gy, gz);
  return g / std::sqrt(gx*gx + gy*gy + gz*gz);
};

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0;
  const FT box_c = (argc > 2) ? std::abs(std::stod(argv[2])) : 1.;
  const std::size_t grid_n = (argc > 3) ? std::stoi(argv[3]) : 50;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox { -box_c, -box_c, -box_c, box_c, box_c, box_c };
  Grid grid { bbox, CGAL::make_array<std::size_t>(grid_n, grid_n, grid_n) };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  // fill up values & gradients
  Values values { devil_value, grid };
  Gradients gradients { devil_gradient, grid };

  // Below is equivalent to:
  //   Domain domain { grid, values, gradients };
  Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;
  CGAL::Isosurfacing::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles,
                                                                       CGAL::parameters::do_not_triangulate_faces(true)
                                                                                        .constrain_to_cell(false));

  std::cout << "Soup #vertices: " << points.size() << std::endl;
  std::cout << "Soup #triangles: " << triangles.size() << std::endl;

  if(!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles)) {
    std::cerr << "Warning: the soup is not a 2-manifold surface, non-manifoldness?..." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Create the mesh..." << std::endl;
  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, mesh);

  CGAL::IO::write_polygon_mesh("dual_contouring.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
