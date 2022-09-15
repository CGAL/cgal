#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {

    Grid grid(100, 100, 100, {-1, -1, -1, 1, 1, 1});

    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                const FT pos_x = x * grid.get_spacing()[0] + grid.get_bbox().xmin();
                const FT pos_y = y * grid.get_spacing()[1] + grid.get_bbox().ymin();
                const FT pos_z = z * grid.get_spacing()[2] + grid.get_bbox().zmin();

                const Vector direction(pos_x, pos_y, pos_z);
                const FT distance = CGAL::approximate_sqrt(direction.squared_length());

                grid.value(x, y, z) = distance;
                grid.gradient(x, y, z) = direction / distance;
            }
        }
    }

    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
