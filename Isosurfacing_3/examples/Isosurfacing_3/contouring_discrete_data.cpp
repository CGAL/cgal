#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_gradients_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
using Gradients = CGAL::Isosurfacing::Interpolated_discrete_gradients_3<Grid>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

void run_marching_cubes(const Grid& grid,
                        const FT isovalue)
{
  using Domain = IS::Marching_cubes_domain_3<Grid, Values, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;

  // fill up values
  Values values { grid };

  for(std::size_t i=0; i<grid.xdim(); ++i) {
    for(std::size_t j=0; j<grid.ydim(); ++j) {
      for(std::size_t k=0; k<grid.zdim(); ++k)
      {
        const Point& p = grid.point(i,j,k);
        const FT d = sqrt(CGAL::squared_distance(p, Point(CGAL::ORIGIN)));
        values(i,j,k) = d;
      }
    }
  }

  Domain domain { grid, values };

  Point_range points;
  Polygon_range triangles;

  // run marching cubes isosurfacing
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("marching_cubes_discrete.off", points, triangles);
}

void run_dual_contouring(const Grid& grid,
                         const FT isovalue)
{
  using Domain = IS::Dual_contouring_domain_3<Grid, Values, Gradients, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  // fill up values and gradients
  Values values { grid };
  Gradients gradients { grid };

  for(std::size_t i=0; i<grid.xdim(); ++i) {
    for(std::size_t j=0; j<grid.ydim(); ++j) {
      for(std::size_t k=0; k<grid.zdim(); ++k)
      {
        const Point& p = grid.point(i,j,k);
        const FT d = sqrt(CGAL::squared_distance(p, Point(CGAL::ORIGIN)));
        values(i,j,k) = d;

        if(d != 0)
          gradients(i,j,k) = Vector(CGAL::ORIGIN, p) / d;
        else
          gradients(i,j,k) = CGAL::NULL_VECTOR;
      }
    }
  }

  Domain domain { grid, values, gradients };

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_discrete.off", points, triangles);
}

int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.8;

  // create bounding box and grid
  const CGAL::Bbox_3 bbox { -1., -1., -1., 1., 1., 1. };
  Grid grid { bbox, CGAL::make_array<std::size_t>(30, 30, 30) };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  run_marching_cubes(grid, isovalue);

  run_dual_contouring(grid, isovalue);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
