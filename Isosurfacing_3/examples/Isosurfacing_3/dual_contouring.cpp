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

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.8;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox { -1., -1., -1., 1., 1., 1. };
  Grid grid { bbox, CGAL::make_array<std::size_t>(30, 30, 30) };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  // fill up values and gradients
  auto sphere_value_fn = [](const Point& p) -> FT
  {
    return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
  };

  auto sphere_gradient_fn = [](const Point& p) -> Vector
  {
    const FT d = sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
    return Vector(CGAL::ORIGIN, p) / d;
  };

  Values values { sphere_value_fn, grid };
  Gradients gradients { sphere_gradient_fn, grid };

  // Below is equivalent to:
  //   Domain domain { grid, values, gradients };
  Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;
  CGAL::Isosurfacing::dual_contouring(domain, isovalue, points, triangles,
                                      CGAL::parameters::do_not_triangulate_faces(true)
                                                       .constrain_to_cell(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring.off", points, triangles);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
