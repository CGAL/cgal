#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Default_gradients.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Explicit_cartesian_grid_domain.h>
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

    // create bounding box and grid
    const CGAL::Bbox_3 bbox(-1.0, -1.0, -1.0,  1.0, 1.0, 1.0);
    std::shared_ptr<Grid> grid = std::make_shar ed<Grid>(30, 30, 30, bbox);

    // compute field values and gradients
    for (std::size_t x = 0; x < grid->xdim(); x++) {
        for (std::size_t y = 0; y < grid->ydim(); y++) {
            for (std::size_t z = 0; z < grid->zdim(); z++) {

                const FT pos_x = x * grid->get_spacing()[0] + bbox.xmin();
                const FT pos_y = y * grid->get_spacing()[1] + bbox.ymin();
                const FT pos_z = z * grid->get_spacing()[2] + bbox.zmin();

                const Vector direction(pos_x, pos_y, pos_z);
                const FT distance = CGAL::approximate_sqrt(direction.squared_length());

                grid->value(x, y, z) = distance;
                grid->gradient(x, y, z) = direction / distance; // TODO: check division / 0
            }
        }
    }

    // gradient field
    CGAL::Isosurfacing::Explicit_cartesian_grid_gradient<Kernel> gradient(grid);

    // create domain from scalar and gradient fields
    auto domain = CGAL::Isosurfacing::create_explicit_cartesian_grid_domain<Kernel>(grid, gradient);

    Point_range points;
    Polygon_range polygons;

    // run dual contouring isosurfacing
    CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

    // write output indexed surface mesh to file, in OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
