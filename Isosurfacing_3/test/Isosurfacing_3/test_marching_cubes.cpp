#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <string>

using K = CGAL::Simple_cartesian<double>;
using FT = typename K::FT;
using Point = typename K::Point_3;
using Vector = typename K::Vector_3;

using Point_range = std::vector<Point>;
using Triangle_range = std::vector<std::array<std::size_t, 3> >;

using Mesh = CGAL::Surface_mesh<Point>;

#define CGAL_TESTUISTE_ISOSURFACING_OUTPUT

namespace IS = CGAL::Isosurfacing;

struct Sphere_function
{
  FT operator()(const Point& p) const
  {
    return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
  }
};

void test_implicit_sphere()
{
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Value_function_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  const CGAL::Bbox_3 bbox {-1., -1., -1., 1., 1., 1.};
  const Vector spacing { 0.1, 0.1, 0.1 };
  const FT isovalue = 0.8;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes (Implicit sphere) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };
  Values values { Sphere_function(), grid };
  Domain domain { grid, values };

  Point_range points;
  Triangle_range triangles;
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << triangles.size() << std::endl;

#ifdef CGAL_TESTUISTE_ISOSURFACING_OUTPUT
  CGAL::IO::write_polygon_soup("MC_implicit_sphere.off", points, triangles);
#endif

  assert(points.size() && triangles.size());
  assert(!has_duplicate_points(points, triangles));
  assert(!has_duplicate_polygons(points, triangles));
  assert(!has_isolated_vertices(points, triangles));

  assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
  Mesh m;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

  assert(is_manifold(m));
  assert(!has_degenerate_faces(m));
}

void test_grid_sphere(const std::size_t n)
{
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  const CGAL::Bbox_3 bbox = {-1., -1., -1., 1., 1., 1.};
  const Vector spacing(2. / (n - 1), 2. / (2 * (n - 1)), 2. / (3*(n - 1)));
  const FT isovalue = 0.777;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes (Implicit sphere) with n = " << n << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  Values values { grid };
  Sphere_function sphere_function;
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const Point pos(x * spacing.x() + bbox.xmin(),
                        y * spacing.y() + bbox.ymin(),
                        z * spacing.z() + bbox.zmin());

        values(x, y, z) = sphere_function(pos);
      }
    }
  }

  Domain domain { grid, values };

  Point_range points;
  Triangle_range triangles;
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << triangles.size() << std::endl;

#ifdef CGAL_TESTUISTE_ISOSURFACING_OUTPUT
  const std::string test_name = "test_grid_sphere(" + std::to_string(n) + ")";
  CGAL::IO::write_polygon_soup(test_name + ".off", points, triangles);
#endif

  assert(!has_duplicate_points(points, triangles));
  assert(!has_duplicate_polygons(points, triangles));
  assert(!has_isolated_vertices(points, triangles));

  assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
  Mesh m;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

  assert(is_manifold(m));
  assert(!has_degenerate_faces(m));
}

int main(int, char**)
{
  test_implicit_sphere();
  test_grid_sphere(2);
  test_grid_sphere(3);
  test_grid_sphere(10);
  test_grid_sphere(50);
  test_grid_sphere(100);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
