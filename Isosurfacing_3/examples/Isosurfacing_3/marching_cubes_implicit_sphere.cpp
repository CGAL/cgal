
#include <CGAL/Isosurfacing_3/Implicit_domain.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    // distance function of a sphere at the origin
    auto sphere_function = [](const Point_3 &point) {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    };

    // create a domain with bounding box [-1, 1]^3 and grid spacing 0.02
    CGAL::Implicit_domain<Kernel, decltype(sphere_function)> domain(sphere_function, CGAL::Bbox_3(-1, -1, -1, 1, 1, 1),
                                                                    Vector_3(0.02f, 0.02f, 0.02f));

    // prepare collections for the resulting polygon soup
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0.8
    CGAL::make_triangle_mesh_using_marching_cubes(domain, 0.8f, points, polygons);

    // save the polygon soup in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}