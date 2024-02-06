#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_gradient_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  // create bounding box and grid
  const CGAL::Bbox_3 bbox{-1., -1., -1.,  1., 1., 1.};
  Grid grid { 30, 30, 30, bbox };

  // compute field values and gradients
  for(std::size_t x = 0; x < grid.xdim(); ++x) {
    for(std::size_t y = 0; y < grid.ydim(); ++y) {
      for(std::size_t z = 0; z < grid.zdim(); ++z)
      {
        const FT pos_x = x * grid.spacing()[0] + bbox.xmin();
        const FT pos_y = y * grid.spacing()[1] + bbox.ymin();
        const FT pos_z = z * grid.spacing()[2] + bbox.zmin();

        const Vector direction(pos_x, pos_y, pos_z);
        const FT distance = CGAL::approximate_sqrt(direction.squared_length());

        grid.value(x, y, z) = distance;

        if(distance != 0)
          grid.gradient(x, y, z) = direction / distance;
        else
          grid.gradient(x, y, z) = CGAL::NULL_VECTOR;
      }
    }
  }

  // gradient field
  CGAL::Isosurfacing::Explicit_Cartesian_grid_gradient_3<Grid> gradient(grid);

  // create domain from scalar and gradient fields
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid, gradient);

  Point_range points;
  Polygon_range polygons;

  // run dual contouring isosurfacing
  CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

  // write output indexed surface mesh to file, in OFF format
  CGAL::IO::write_OFF("dual_contouring_Cartesian_grid.off", points, polygons);

  return EXIT_SUCCESS;
}
