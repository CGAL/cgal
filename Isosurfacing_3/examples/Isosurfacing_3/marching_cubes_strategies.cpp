#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;

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

  // fill up values
  auto sphere_value_fn = [](const Point& p) -> FT
  {
    return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
  };

  Values values { sphere_value_fn, grid };
  Domain domain { grid, values };
  // the domain could also be created with:
  // auto domain = CGAL::Isosurfacing::create_marching_cubes_domain_3(grid, values);

  // MC base version
  {
    Point_range points;
    Polygon_range triangles;

    CGAL::Real_timer timer;
    timer.start();

    // run marching cubes isosurfacing
    std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;
    CGAL::Isosurfacing::marching_cubes(domain, isovalue, points, triangles,
                                      CGAL::parameters::use_topologically_correct_marching_cubes(false));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("marching_cubes.off", points, triangles);
  }

  // MC topologically correct version
  {
    Point_range points;
    Polygon_range triangles;

    CGAL::Real_timer timer;
    timer.start();

    // run marching cubes isosurfacing
    std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;
    CGAL::Isosurfacing::marching_cubes(domain, isovalue, points, triangles,
                                      CGAL::parameters::use_topologically_correct_marching_cubes(true));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("marching_cubes_TMC.off", points, triangles);
  }

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
