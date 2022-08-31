#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {

    const std::string fname = "../data/skull_2.9.inr";
    // Load image
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }

    Grid grid(image);

    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, 2.9f, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
