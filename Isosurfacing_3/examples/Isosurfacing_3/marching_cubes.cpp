#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <cmath>
#include <iostream>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

using Mesh = CGAL::Surface_mesh<Point>;

auto value_fn = [](const Point& p) -> FT
{
  const FT& x = p.x(), y = p.y(), z = p.z();

  // "Klein Bottle" - https://www-sop.inria.fr/galaad/surface/
  return -1-4*y*z*z*x*x-2*y+6*z*z*x*x*y*y-16*x*z+16*x*z*y*y+3*x*x+7*y*y+11*z*z-11*z*z*z*z+z*z*z*z*z*z-3*x*x*x*x-7*y*y*y*y+x*x*x*x*x*x+y*y*y*y*y*y-14*z*z*x*x-18*z*z*y*y+3*z*z*z*z*x*x+3*z*z*z*z*y*y-10*x*x*y*y-4*y*y*y*z*z+3*z*z*x*x*x*x+3*z*z*y*y*y*y+16*x*x*x*z+3*x*x*x*x*y*y+3*x*x*y*y*y*y+4*x*x*y-12*z*z*y-2*x*x*x*x*y-4*x*x*y*y*y-2*z*z*z*z*y+16*x*z*z*z+12*y*y*y-2*y*y*y*y*y-32*x*z*y;
};

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.1;
  const FT box_c = (argc > 2) ? std::abs(std::stod(argv[2])) : 5.;
  const std::size_t grid_n = (argc > 3) ? std::stoi(argv[3]) : 50;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox { -box_c, -box_c, -box_c, box_c, box_c, box_c };
  Grid grid { bbox, CGAL::make_array<std::size_t>(grid_n, grid_n, grid_n) };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  // fill up values
  Values values { value_fn, grid };

  // Below is equivalent to:
  //   Domain domain { grid, values };
  Domain domain = CGAL::Isosurfacing::create_marching_cubes_domain_3(grid, values);

  Point_range points;
  Polygon_range triangles;

  // run marching cubes isosurfacing
  std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;
  CGAL::Isosurfacing::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles,
                                                                      CGAL::parameters::use_topologically_correct_marching_cubes(true));

  std::cout << "Soup #vertices: " << points.size() << std::endl;
  std::cout << "Soup #triangles: " << triangles.size() << std::endl;

  if(!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles)) {
    std::cerr << "Warning: the soup is not a 2-manifold surface, non-manifoldness?..." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Convert the soup to a triangle mesh..." << std::endl;
  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, mesh);

  CGAL::IO::write_polygon_mesh("marching_cubes.off", mesh, CGAL::parameters::stream_precision(17));

  // Let's remesh it to something nicer looking
  std::cout << "Remeshing..." << std::endl;
  const FT target_edge_length = box_c / 50;
  CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
                                                     CGAL::parameters::number_of_iterations(5)
                                                                      .number_of_relaxation_steps(5));

  CGAL::IO::write_polygon_mesh("marching_cubes-remeshed.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
