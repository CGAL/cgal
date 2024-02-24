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

#include <CGAL/IO/polygon_soup_io.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

auto implicit_function = [](const Point& q) -> FT
{
  auto cyl = [](const Point& q) { return IS::Shapes::infinite_cylinder<Kernel>(Point(0,0,0), Vector(0,0,1), 0.5, q); };
  auto cube = [](const Point& q) { return IS::Shapes::box<Kernel>(Point(-0.5,-0.5,-0.5), Point(0.5,0.5,0.5), q); };
  auto cyl_and_cube = [&](const Point& q) { return IS::Shapes::shape_union<Kernel>(cyl, cube, q); };

  auto sphere = [](const Point& q) { return IS::Shapes::sphere<Kernel>(Point(0,0,0.5), 0.75, q); };
  return IS::Shapes::shape_difference<Kernel>(cyl_and_cube, sphere, q);
};

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.;
  std::cout << "Isovalue: " << isovalue << std::endl;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox = {-2., -2., -2., 2., 2., 2.};
  Grid grid { bbox, CGAL::make_array<std::size_t>(150, 150, 150) };

  std::cout << "Bbox: " << grid.bbox() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  // fill up values and gradients
  Values values { implicit_function, grid };
  Gradients gradients { values };

  const bool triangulate_faces = false;

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
    CGAL::IO::write_polygon_soup("dual_contouring_implicit_parallel.off", points, triangles);
  }

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
    CGAL::IO::write_polygon_soup("dual_contouring_implicit_sequential.off", points, triangles);
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
    CGAL::IO::write_polygon_soup("marching_cubes_implicit_sequential.off", points, triangles);
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
    CGAL::IO::write_polygon_soup("marching_cubes_implicit_parallel.off", points, triangles);
  }

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
