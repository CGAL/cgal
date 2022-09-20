
#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_domains.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    // create a cartesian grid with 100^3 grid points and the bounding box [-1, 1]^3
    Grid grid(50, 50, 50, {-1, -1, -1, 1, 1, 1});

    // calculate the value at all grid points
    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                const FT pos_x = x * grid.get_spacing()[0] + grid.get_bbox().xmin();
                const FT pos_y = y * grid.get_spacing()[1] + grid.get_bbox().ymin();
                const FT pos_z = z * grid.get_spacing()[2] + grid.get_bbox().zmin();

                // distance to the origin
                grid.value(x, y, z) = std::sqrt(pos_x * pos_x + pos_y * pos_y + pos_z * pos_z);
            }
        }
    }

    // create a domain from the grid
    auto domain = CGAL::Isosurfacing::create_explicit_cartesian_grid_domain<Kernel>(grid);

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0.8
    CGAL::Isosurfacing::marching_cubes(domain, 0.8f, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
