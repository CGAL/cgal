#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_domains.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {

    const std::string fname = CGAL::data_file_path("skull_2.9.inr");  // TODO: get other examples with rights

    // load the image
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }

    // convert it to a cartesian grid
    const Grid grid(image);

    // create a domain from the grid
    auto domain = CGAL::Isosurfacing::create_explicit_cartesian_grid_domain<Kernel>(grid);

    // prepare collections for the result
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 2.9
    CGAL::Isosurfacing::marching_cubes(domain, 2.9, points, polygons);

    // save the result in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
