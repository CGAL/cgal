#include <CGAL/Simple_cartesian.h>
#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  // 3D bounding box [-1, 1]^3 and Cartesian grid
  const CGAL::Bbox_3 bbox{-1., -1., -1.,  1., 1., 1.};
  Grid grid { 50, 50, 50, bbox };

  // compute and store function values at all grid points
  for(std::size_t x = 0; x < grid.xdim(); ++x) {
    for(std::size_t y = 0; y < grid.ydim(); ++y) {
      for(std::size_t z = 0; z < grid.zdim(); ++z)
      {
        const FT pos_x = x * grid.spacing()[0] + bbox.xmin();
        const FT pos_y = y * grid.spacing()[1] + bbox.ymin();
        const FT pos_z = z * grid.spacing()[2] + bbox.zmin();

        // Euclidean distance to the origin
        grid.value(x, y, z) = sqrt(pos_x * pos_x + pos_y * pos_y + pos_z * pos_z);
      }
    }
  }

  // create a domain from the grid
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // prepare collections for the result
  Point_range points;
  Polygon_range triangles;

  // run marching cubes with an isovalue of 0.8
  CGAL::Isosurfacing::marching_cubes(domain, 0.8, points, triangles);

  // save output indexed surface mesh to file, in the OFF format
  CGAL::IO::write_OFF("marching_cubes_Cartesian_grid_sphere.off", points, triangles);

  return EXIT_SUCCESS;
}
