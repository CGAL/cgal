#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <iostream>

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef CGAL::Surface_mesh<Point> Mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;


inline Kernel::FT distance_to_mesh(const Tree& tree, const Point& p) {
    const Point& x = tree.closest_point(p);
    return std::sqrt((p - x).squared_length());
}

int main() {
    const std::string input_name = "../../../../data/bunny.off";
    const int n_voxels = 50;
    const FT offset_value = 0.02;

    Mesh mesh_input;
    if (!CGAL::IO::read_OFF(input_name, mesh_input)) {
        std::cerr << "Could not read mesh" << std::endl;
        exit(-1);
    }

    // compute bounding box
    CGAL::Bbox_3 aabb_grid = CGAL::Polygon_mesh_processing::bbox(mesh_input);
    Vector aabb_increase_vec = Vector(offset_value + 0.01, offset_value + 0.01, offset_value + 0.01);
    aabb_grid += (Point(aabb_grid.xmax(), aabb_grid.ymax(), aabb_grid.zmax()) + aabb_increase_vec).bbox();
    aabb_grid += (Point(aabb_grid.xmin(), aabb_grid.ymin(), aabb_grid.zmin()) - aabb_increase_vec).bbox();

    // construct AABB tree
    Tree tree(mesh_input.faces_begin(), mesh_input.faces_end(), mesh_input);

    CGAL::Side_of_triangle_mesh<Mesh, CGAL::GetGeomTraits<Mesh>::type> sotm(mesh_input);

    Grid grid(n_voxels, n_voxels, n_voxels, aabb_grid);

    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    for (std::size_t z = 0; z < grid.zdim(); z++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t x = 0; x < grid.xdim(); x++) {

                const FT pos_x = x * grid.voxel_x() + grid.offset_x();
                const FT pos_y = y * grid.voxel_y() + grid.offset_y();
                const FT pos_z = z * grid.voxel_z() + grid.offset_z();
                const Point p(pos_x, pos_y, pos_z);

                grid.value(x, y, z) = distance_to_mesh(tree, p);

                const bool is_inside = (sotm(p) == CGAL::ON_BOUNDED_SIDE);
                if (is_inside) {
                    grid.value(x, y, z) *= -1;
                }
            }
        }
    }

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, offset_value, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
