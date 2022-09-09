
#include <CGAL/Implicit_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <tbb/concurrent_vector.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef tbb::concurrent_vector<Point> Point_range;
typedef tbb::concurrent_vector<std::vector<std::size_t>> Polygon_range;

int main() {
    // distance to the origin
    auto sphere_function = [](const Point& point) {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    };

    // create a domain with bounding box [-1, 1]^3 and grid spacing 0.02
    CGAL::Isosurfacing::Implicit_domain<Kernel, decltype(sphere_function), decltype(CGAL::Isosurfacing::Default_gradient<Kernel, decltype(sphere_function)>(sphere_function))> domain(
        {-1, -1, -1, 1, 1, 1}, Vector(0.02f, 0.02f, 0.02f), sphere_function, CGAL::Isosurfacing::Default_gradient<Kernel, decltype(sphere_function)>(sphere_function));  // TODO: this is ugly

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0.8
    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, 0.8f, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
