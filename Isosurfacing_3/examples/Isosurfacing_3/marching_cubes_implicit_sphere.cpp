#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Implicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>

#include <CGAL/Bbox_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <vector>

#include <tbb/tick_count.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Vector = typename Kernel::Vector_3;
using Point = typename Kernel::Point_3;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  const CGAL::Bbox_3 bbox { -1.0, -1.0, -1.0,  1.0, 1.0, 1.0 };
  const FT spacing = 0.002;
  const Vector vec_spacing(spacing, spacing, spacing);

  // Euclidean distance function to the origin
  auto sphere_function = [&](const Point& p) -> FT
  {
    return sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
  };

  // create a domain with given bounding box and grid spacing
  auto domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain<Kernel>(bbox, vec_spacing, sphere_function);

  // prepare collections for the output indexed mesh
  Point_range points;
  Polygon_range polygons;

  const tbb::tick_count start = tbb::tick_count::now();

  // execute marching cubes with an isovalue of 0.8
  CGAL::Isosurfacing::marching_cubes(domain, 0.8, points, polygons);

  const tbb::tick_count end = tbb::tick_count::now();
  std::cout << (end - start).seconds() << std::endl;

  // save ouput indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
