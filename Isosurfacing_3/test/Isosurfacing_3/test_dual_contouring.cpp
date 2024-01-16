#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/Implicit_Cartesian_grid_domain_3.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include "Timer.h"

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;
using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  const Vector spacing(0.002, 0.002, 0.02);
  const CGAL::Bbox_3 bbox = {-1, -1, -1, 1, 1, 1};

  auto sphere_function = [](const Point& point) -> FT
  {
    return sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
  };

  using Gradient = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel, decltype(sphere_function)>;

  auto implicit_domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain<Kernel>(
      bbox, spacing, sphere_function, Gradient(sphere_function, 0.0001));

  const std::size_t nx = static_cast<std::size_t>(2.0 / spacing.x());
  const std::size_t ny = static_cast<std::size_t>(2.0 / spacing.y());
  const std::size_t nz = static_cast<std::size_t>(2.0 / spacing.z());

  Grid grid { nx, ny, nz, bbox };

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const Point pos(x * spacing.x() + bbox.xmin(),
                        y * spacing.y() + bbox.ymin(),
                        z * spacing.z() + bbox.zmin());

        grid.value(x, y, z) = sphere_function(pos);
      }
    }
  }

  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/skull_2.9.inr");

  // Load image
  // CGAL::Image_3 image;
  // if(!image.read(fname)) {
  //    std::cerr << "Error: Cannot read file " << fname << std::endl;
  //    return EXIT_FAILURE;
  //}
  // Grid grid(image);

  auto grid_domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  Point_range points;
  Polygon_range polygons;

  {
    ScopeTimer timer;
    CGAL::Isosurfacing::dual_contouring<CGAL::Parallel_tag>(implicit_domain, 0.8f, points, polygons);
  }

  // @todo compare results with mesh_3

  // Mesh mesh;
  // CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  // CGAL::IO::write_OFF("result.off", mesh);
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
