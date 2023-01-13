#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Implicit_Cartesian_grid_domain_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/IO/OFF.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Vector = typename Kernel::Vector_3;
using Point = typename Kernel::Point_3;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  const FT alpha = 5.01;

  auto iwp_value = [alpha](const Point& point)
  {
    const FT x = alpha * (point.x() + 1) * CGAL_PI;
    const FT y = alpha * (point.y() + 1) * CGAL_PI;
    const FT z = alpha * (point.z() + 1) * CGAL_PI;
    return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z);  // isovalue = 0
  };

  auto iwp_gradient = [alpha](const Point& point)
  {
    const FT x = alpha * (point.x() + 1) * CGAL_PI;
    const FT y = alpha * (point.y() + 1) * CGAL_PI;
    const FT z = alpha * (point.z() + 1) * CGAL_PI;

    const FT gx = CGAL_PI * alpha * sin(x) * (cos(y) * (cos(z) - 1.0) - cos(z));
    const FT gy = CGAL_PI * alpha * sin(y) * (cos(x) * (cos(z) - 1.0) - cos(z));
    const FT gz = CGAL_PI * alpha * sin(z) * (cos(x) * (cos(y) - 1.0) - cos(y));
    return Vector(gx, gy, gz);
  };

  const CGAL::Bbox_3 bbox{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
  const FT spacing = 0.5;
  const Vector vec_spacing(spacing, spacing, spacing);

  // create a domain with given bounding box and grid spacing
  auto domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain<Kernel>(bbox, vec_spacing, iwp_value, iwp_gradient);

  // prepare collections for the result
  Point_range points;
  Polygon_range polygons;

  // run marching cubes with an isovalue of 0.0
  CGAL::Isosurfacing::dual_contouring(domain, 0.0, points, polygons);

  // save the result in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
