#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <tbb/concurrent_vector.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef tbb::concurrent_vector<Point> Point_range;
typedef tbb::concurrent_vector<std::vector<std::size_t>> Polygon_range;

int main() {

    const std::string fname = "../../../data/skull_2.9.inr";

    // load the image
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }

    // convert it to a cartesian grid
    const Grid grid(image);

    // create a domain from the grid
    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 2.9
    CGAL::Isosurfacing::make_quad_mesh_using_dual_contouring(domain, 2.9, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
