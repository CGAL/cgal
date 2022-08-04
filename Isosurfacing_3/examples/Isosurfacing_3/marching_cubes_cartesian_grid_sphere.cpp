
#include <CGAL/Isosurfacing_3/Cartesian_grid_domain.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    // create a cartesian grid with 100^3 grid points and the bounding box [-1, 1]^3
    CGAL::Cartesian_grid_3<Kernel> grid(100, 100, 100, CGAL::Bbox_3(-1, -1, -1, 1, 1, 1));

    // calculate the value at all grid points
    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                // distance function of a sphere at the origin
                grid.value(x, y, z) = std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z())
            }
        }
    }

    // create a domain with bounding box [-1, 1]^3 and grid spacing 0.02
    CGAL::Cartesian_grid_domain<Kernel> domain(grid);

    // prepare collections for the resulting polygon soup
    Point_range points;
    Polygon_range polygons;

    // execute marching cubes with an isovalue of 0.8
    CGAL::make_triangle_mesh_using_marching_cubes(domain, 0.8f, points, polygons);

    // convert the polygon soup to a surface mesh
    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

    // save the mesh in the OFF format
    CGAL::IO::write_OFF("result.off", mesh);
}