#include <CGAL/Bbox_3.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Implicit_cartesian_grid_domain.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

int main() {
    const FT alpha = 5.01;

    auto iwp_value = [alpha](const Point& point) {
        const FT x = alpha * (point.x() + 1) * M_PI;
        const FT y = alpha * (point.y() + 1) * M_PI;
        const FT z = alpha * (point.z() + 1) * M_PI;
        return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z);  // iso-value = 0
    };

    auto iwp_gradient = [alpha](const Point& point) {
        const FT x = alpha * (point.x() + 1) * M_PI;
        const FT y = alpha * (point.y() + 1) * M_PI;
        const FT z = alpha * (point.z() + 1) * M_PI;

        const FT gx = M_PI * alpha * sin(x) * (cos(y) * (cos(z) - 1.0) - cos(z));
        const FT gy = M_PI * alpha * sin(y) * (cos(x) * (cos(z) - 1.0) - cos(z));
        const FT gz = M_PI * alpha * sin(z) * (cos(x) * (cos(y) - 1.0) - cos(y));
        return Vector(gx, gy, gz);
    };

    const CGAL::Bbox_3 bbox{-1, -1, -1, 1, 1, 1};
    const Vector spacing(0.05f, 0.05f, 0.05f);

    // create a domain with bounding box [-1, 1]^3 and grid spacing 0.02
    auto domain =
        CGAL::Isosurfacing::create_implicit_cartesian_grid_domain<Kernel>(bbox, spacing, iwp_value, iwp_gradient);

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0
    CGAL::Isosurfacing::dual_contouring(domain, 0.0f, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
