#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>

#include <CGAL/IO/polygon_soup_io.h>

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

// https://www-sop.inria.fr/galaad/surface/
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
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox { -1, -1, -1, 1, 1, 1 };
  Grid grid { bbox, CGAL::make_array<std::size_t>(50, 50, 50) };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  Values values { devil_value, grid };
  Gradients gradients { devil_gradient, grid };

  // Below is equivalent to:
  //   Domain domain { grid, values, gradients };
  Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;
  CGAL::Isosurfacing::dual_contouring(domain, isovalue, points, triangles,
                                      CGAL::parameters::do_not_triangulate_faces(true)
                                                       .constrain_to_cell(false));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring.off", points, triangles);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
