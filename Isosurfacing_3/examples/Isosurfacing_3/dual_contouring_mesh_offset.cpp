#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Implicit_cartesian_grid_domain.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;

typedef CGAL::Surface_mesh<Point> Mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t> > Polygon_range;

int main() {
    const std::string input_name = CGAL::data_file_path("meshes/cross.off");
    const Vector grid_spacing(0.1, 0.1, 0.1);
    const FT offset_value = 0.2;

    Mesh mesh_input;
    if (!CGAL::IO::read_OFF(input_name, mesh_input)) {
        std::cerr << "Could not read input mesh" << std::endl;
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

    // functors for addressing distance and normal queries
    auto mesh_distance = [&tree](const Point& p) {
        const Point& x = tree.closest_point(p);
        return std::sqrt((p - x).squared_length());
    };

    auto mesh_normal = [&tree](const Point& p) {
        const Point& x = tree.closest_point(p);
        const Vector n = p - x; // TODO: address case where norm = zero
        return n / std::sqrt(n.squared_length()); // normalize output vector
    };

    // create a domain with given bounding box and grid spacing
    auto domain = CGAL::Isosurfacing::create_implicit_cartesian_grid_domain<Kernel>(aabb_grid, grid_spacing,
                                                                                    mesh_distance, mesh_normal);
    // containers for output indexed surface mesh
    Point_range points;
    Polygon_range polygons;

    // dual contouring
    CGAL::Isosurfacing::dual_contouring(domain, offset_value, points, polygons);

    // save output indexed mesh to a file, in the OFF format
    CGAL::IO::write_OFF("result.off", points, polygons);
}
