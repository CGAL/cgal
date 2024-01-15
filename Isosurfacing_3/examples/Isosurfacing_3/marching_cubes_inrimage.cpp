#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = typename Kernel::Point_3;
using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/skull_2.9.inr");

  // load volumetric image from a file
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read image file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  Grid grid{image};

  for (std::size_t i = 0; i < grid.xdim(); i++)
    for (std::size_t j = 0; j < grid.ydim(); j++)
      for (std::size_t k = 0; k < grid.zdim(); k++)
        grid.value(i, j, k) = 2 * 1120 - grid.value(i, j, k);

  // create a domain from the grid
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // prepare collections for the output indexed soup
  Point_range points;
  Polygon_range polygons;

  // execute marching cubes
  CGAL::Isosurfacing::marching_cubes(domain, 1120 /*isovalue*/, points, polygons);

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("marching_cubes_inrimage.off", points, polygons);

  return EXIT_SUCCESS;
}
