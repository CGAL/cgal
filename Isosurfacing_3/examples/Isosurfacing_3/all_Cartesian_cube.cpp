#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

// return 1.0 if value has positive sign, and -1.0 otherwise
FT sign(FT value)
{
  return (value > 0.0) - (value < 0.0);
}

int main(int, char**)
{
  // create a Cartesian grid with 7^3 grid points and the bounding box [-1, 1]^3
  const CGAL::Bbox_3 bbox{-1., -1., -1.,  1., 1., 1.};
  Grid grid { 7, 7, 7, bbox };

  // calculate the value at all grid points
  for(std::size_t i=0; i<grid.xdim(); ++i) {
    for(std::size_t j=0; j<grid.ydim(); ++j) {
      for(std::size_t k=0; k<grid.zdim(); ++k)
      {
        const FT pos_x = i * grid.spacing()[0] + bbox.xmin();
        const FT pos_y = j * grid.spacing()[1] + bbox.ymin();
        const FT pos_z = k * grid.spacing()[2] + bbox.zmin();

        // L_inf distance to the origin
        grid.value(i, j, k) = (std::max)({std::abs(pos_x), std::abs(pos_y), std::abs(pos_z)});
      }
    }
  }

  // compute function gradient
  auto cube_gradient = [](const Point& p)
  {
    // the normal depends on the side of the cube
    const FT max_value = std::max({ std::abs(p.x()), std::abs(p.y()), std::abs(p.z())});

    Vector g(0.0, 0.0, 0.0);
    if(max_value == std::abs(p.x()))
      g += Vector(sign(p.x()), 0.0, 0.0);

    if(max_value == std::abs(p.y()))
      g += Vector(0.0, sign(p.y()), 0.0);

    if(max_value == std::abs(p.z()))
      g += Vector(0.0, 0.0, sign(p.z()));

    const FT length_sq = g.squared_length();
    if(length_sq > 0.00001)
      g /= CGAL::approximate_sqrt(length_sq);

    return g;
  };

  // create domain from given grid and gradient
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid, cube_gradient);

  // containers for output indexed surface meshes
  Point_range points_mc, points_dc;
  Polygon_range polygons_mc, polygons_dc;

  // run topologically correct Marching Cubes and Dual Contouring with given isovalue
  const FT isovalue = 0.88;
  CGAL::Isosurfacing::marching_cubes(domain, isovalue, points_mc, polygons_mc);
  CGAL::Isosurfacing::dual_contouring(domain, isovalue, points_dc, polygons_dc);

  // save output indexed meshes to files, in the OFF format @fixme these are not meshes...
  CGAL::IO::write_OFF("output_mc.off", points_mc, polygons_mc);
  CGAL::IO::write_OFF("output_dc.off", points_dc, polygons_dc);

  return EXIT_SUCCESS;
}
