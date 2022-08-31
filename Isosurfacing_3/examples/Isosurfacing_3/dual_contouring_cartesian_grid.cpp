#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point_3;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {

    Grid grid(100, 100, 100, {-1, -1, -1, 1, 1, 1});

    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                const FT pos_x = x * grid.voxel_x() + grid.offset_x();
                const FT pos_y = y * grid.voxel_y() + grid.offset_y();
                const FT pos_z = z * grid.voxel_z() + grid.offset_z();

                grid.value(x, y, z) = std::sqrt(pos_x * pos_x + pos_y * pos_y + pos_z * pos_z);
            }
        }
    }

    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_quad_mesh_using_dual_contouring(domain, 0.8, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
