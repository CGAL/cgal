
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>

#include "test_util.h"

#include <iostream>
#include <string>

#define WRITE_OFF

struct Sphere_function
{
  FT operator()(const Point& point) const
  {
    return sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
  }
};

template <typename Domain_>
void run(const Domain_& domain,
         const FT isovalue,
         Point_range& points,
         Polygon_range& polygons)
{
  CGAL::Isosurfacing::marching_cubes<CGAL::Parallel_tag>(domain, isovalue, points, polygons);
}

void test_implicit_sphere()
{
  const std::string test_name = "test_implicit_sphere()";

  const Vector spacing(0.2, 0.2, 0.2);
  const CGAL::Bbox_3 bbox = {-1, -1, -1, 1, 1, 1};

  auto domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain<Kernel>(bbox, spacing, Sphere_function());

  Point_range points;
  Polygon_range polygons;
  run(domain, 0.8, points, polygons);

#ifdef WRITE_OFF
  CGAL::IO::write_OFF(test_name + ".off", points, polygons);
#endif

  assert(is_polygon_mesh(polygons));
  Mesh m = to_mesh(points, polygons);

  assert(is_manifold(m));
  assert(!has_degenerate_faces(m));

  std::cout << "Test passed: " << test_name << std::endl;
}

void test_grid_sphere(const std::size_t n)
{
  const std::string test_name = "test_grid_sphere(" + std::to_string(n) + ")";

  const CGAL::Bbox_3 bbox = {-1, -1, -1, 1, 1, 1};
  const Vector spacing(2.0 / (n - 1), 2.0 / (n - 1), 2.0 / (n - 1));

  Sphere_function sphere_function;

  Grid grid { bbox, CGAL::make_array<std::size_t>(n, n, n) };

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

  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  Point_range points;
  Polygon_range polygons;
  run(domain, 0.777, points, polygons);

#ifdef WRITE_OFF
  CGAL::IO::write_OFF(test_name + ".off", points, polygons);
#endif

  assert(is_polygon_mesh(polygons));
  Mesh m = to_mesh(points, polygons);

  assert(is_manifold(m));
  assert(!has_degenerate_faces(m));

  std::cout << "Test passed: " << test_name << std::endl;
}

int main(int, char**)
{
  test_implicit_sphere();
  test_grid_sphere(2);
  test_grid_sphere(3);
  test_grid_sphere(10);
  test_grid_sphere(11);
  test_grid_sphere(100);

  std::cout << "All tests passed" << std::endl;

  return EXIT_SUCCESS;
}
