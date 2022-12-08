#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Explicit_cartesian_grid_domain.h>
#include <CGAL/Implicit_cartesian_grid_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Finite_difference_gradient.h>

#include "Timer.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    const Vector spacing(0.002, 0.002, 0.02);
    const CGAL::Bbox_3 bbox = {-1, -1, -1, 1, 1, 1};

    auto sphere_function = [](const Point& point) {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    };

    typedef CGAL::Isosurfacing::Finite_difference_gradient<Kernel, decltype(sphere_function)> Gradient;

    auto implicit_domain = CGAL::Isosurfacing::create_implicit_cartesian_grid_domain<Kernel>(
        bbox, spacing, sphere_function, Gradient(sphere_function, 0.0001));

    const std::size_t nx = static_cast<std::size_t>(2.0 / spacing.x());
    const std::size_t ny = static_cast<std::size_t>(2.0 / spacing.y());
    const std::size_t nz = static_cast<std::size_t>(2.0 / spacing.z());

    std::shared_ptr<Grid> grid = std::make_shared<Grid>(nx, ny, nz, bbox);

    for (std::size_t x = 0; x < grid->xdim(); x++) {
        for (std::size_t y = 0; y < grid->ydim(); y++) {
            for (std::size_t z = 0; z < grid->zdim(); z++) {

                const Point pos(x * spacing.x() + bbox.xmin(), y * spacing.y() + bbox.ymin(),
                                z * spacing.z() + bbox.zmin());

                grid->value(x, y, z) = sphere_function(pos);
            }
        }
    }

    const std::string fname = "../data/skull_2.9.inr";
    // Load image
    // CGAL::Image_3 image;
    // if (!image.read(fname)) {
    //    std::cerr << "Error: Cannot read file " << fname << std::endl;
    //    return EXIT_FAILURE;
    //}
    // Grid grid(image);

    auto grid_domain = CGAL::Isosurfacing::create_explicit_cartesian_grid_domain<Kernel>(grid);

    Point_range points;
    Polygon_range polygons;

    {
        ScopeTimer timer;
        CGAL::Isosurfacing::dual_contouring<CGAL::Parallel_tag>(implicit_domain, 0.8f, points, polygons);
    }

    // TODO: compare results with mesh_3

    // Mesh mesh;
    // CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

    // CGAL::IO::write_OFF("result.off", mesh);
    CGAL::IO::write_OFF("result.off", points, polygons);

    return 0;
}
