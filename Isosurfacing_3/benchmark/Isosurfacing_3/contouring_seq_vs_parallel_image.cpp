#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/internal/implicit_shapes_helper.h>

#include <CGAL/Real_timer.h>

#include <CGAL/Isosurfacing_3/IO/Image_3.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

int main(int argc, char** argv)
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/skull_2.9.inr");
  const FT isovalue = (argc > 2) ? std::stod(argv[2]) : - 2.9;

  // load volumetric image from a file
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read image file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  Grid grid;
  Values values { grid };
  if(!IS::IO::read_Image_3(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    return EXIT_FAILURE;
  }

  for (std::size_t i=0; i<grid.xdim(); ++i)
    for (std::size_t j=0; j<grid.ydim(); ++j)
      for (std::size_t k=0; k<grid.zdim(); ++k)
        values(i, j, k) = - values(i, j, k); // inside out

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  Gradients gradients { values };

  const bool triangulate_faces = false;

  // DC sequential
  {
    using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;
    Domain domain { grid, values, gradients };

    CGAL::Real_timer timer;
    timer.start();

    Point_range points;
    Polygon_range triangles;

    std::cout << "--- Dual Contouring (Sequential)" << std::endl;
    IS::dual_contouring<CGAL::Sequential_tag>(
      domain, isovalue, points, triangles, CGAL::parameters::do_not_triangulate_faces(!triangulate_faces));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("dual_contouring_image_sequential.off", points, triangles);
  }

  // DC parallel
  {
    using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;
    Domain domain { grid, values, gradients };

    CGAL::Real_timer timer;
    timer.start();

    Point_range points;
    Polygon_range triangles;

    std::cout << "--- Dual Contouring (Parallel)" << std::endl;
    IS::dual_contouring<CGAL::Parallel_if_available_tag>(
      domain, isovalue, points, triangles, CGAL::parameters::do_not_triangulate_faces(!triangulate_faces));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("dual_contouring_image_parallel.off", points, triangles);
  }

  // MC Sequential
  {
    using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;
    Domain domain { grid, values };

    CGAL::Real_timer timer;
    timer.start();

    Point_range points;
    Polygon_range triangles;

    std::cout << "--- Marching Cubes (Sequential)" << std::endl;
    IS::marching_cubes<CGAL::Sequential_tag>(
      domain, isovalue, points, triangles, CGAL::parameters::do_not_triangulate_faces(!triangulate_faces));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("marching_cubes_image_sequential.off", points, triangles);
  }

  // MC parallel
  {
    using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;
    Domain domain { grid, values };

    CGAL::Real_timer timer;
    timer.start();

    Point_range points;
    Polygon_range triangles;

    std::cout << "--- Marching Cubes (Parallel)" << std::endl;
    IS::marching_cubes<CGAL::Parallel_if_available_tag>(
      domain, isovalue, points, triangles, CGAL::parameters::do_not_triangulate_faces(!triangulate_faces));

    timer.stop();

    std::cout << "Output #vertices: " << points.size() << std::endl;
    std::cout << "Output #triangles: " << triangles.size() << std::endl;
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    CGAL::IO::write_polygon_soup("marching_cubes_image_parallel.off", points, triangles);
  }

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
