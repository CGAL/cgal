
#include <CGAL/Bbox_3.h>
#include <CGAL/Implicit_cartesian_grid_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    const CGAL::Bbox_3 bbox{-1, -1, -1, 1, 1, 1};
    const Vector spacing(0.04f, 0.04f, 0.04f);

    auto sphere_function = [&](const Point& p) { return std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()); };

    // create a domain with bounding box [-1, 1]^3 and grid spacing 0.02
    auto domain =
        CGAL::Isosurfacing::create_implicit_cartesian_grid_domain<Kernel>(bbox, spacing, std::move(sphere_function));

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0.8
    CGAL::Isosurfacing::marching_cubes(domain, 0.8f, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
