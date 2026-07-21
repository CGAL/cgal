
#include <Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>

class Mesh_wrapper : public Mesh_smoothing_3::helper_structures::Mixed_element_mesh<GEO::index_t, GEO::index_t, GEO::vec3, GEO::index_range> {
public:
    std::size_t nb_vertices() const override { return mesh.vertices.nb(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const override {
        return mesh.vertices.point(vertex);
    }
    void set_new_vertex_coordinates(GEO::index_t vertex, Point_3 coord) override { mesh.vertices.point(vertex) = coord; }

    GEO::index_range input_cell_range() const override { return GEO::index_range{mesh.cells.begin(), mesh.cells.end()}; }

    Vertex_descriptor get_cell_vertex(GEO::index_t cell, unsigned local_Vertex_descriptor) const override {
        return mesh.cells.vertex(cell, local_Vertex_descriptor);
    };

    GEO::vec3 get_ref_vertex_coordinates(GEO::index_t vertex) const override {
        assert(ref_mesh != nullptr);
        return ref_mesh->vertices.point(vertex);
    }

    Shape const * get_element_shape(GEO::index_t cell) const override {
        switch (mesh.cells.type(cell)) {
            case GEO::MESH_TET:
                return &tet_ref;
                break;
            case GEO::MESH_HEX:
                return &hex_ref;
                break;
            case GEO::MESH_PYRAMID:
                return &py_ref;
                break;
            case GEO::MESH_PRISM:
                return &we_ref;
                break;
            default:
                return nullptr;
                break;
        }
    }

public:
    Mesh_wrapper(GEO::Mesh &mesh_, GEO::Mesh const * ref_mesh_)
    : mesh(mesh_)
    , ref_mesh(ref_mesh_)
    {
        this->has_reference_mesh = (ref_mesh != nullptr);
        we_ref.inverse = true;
        this->assemble(); // CRITICAL
    }

    void set_orientation(bool inv_tet, bool inv_hex, bool inv_pyr, bool inv_wed) {
        tet_ref.inverse = inv_tet;
        hex_ref.inverse = inv_hex;
        py_ref.inverse  = inv_pyr;
        we_ref.inverse  = !inv_wed;
    }

    GEO::Mesh &mesh;
    GEO::Mesh const * ref_mesh;

    Mesh_smoothing_3::Shapes::VTK_TETRAHEDRON<GEO::vec3> tet_ref;
    Mesh_smoothing_3::Shapes::GEOGRAM_HEXAHEDRON<GEO::vec3> hex_ref;
    Mesh_smoothing_3::Shapes::VTK_PYRAMID<GEO::vec3> py_ref;
    Mesh_smoothing_3::Shapes::VTK_WEDGE<GEO::vec3> we_ref;
};

class Boundary_wrapper {
public:
    using Face_descriptor = GEO::index_t;
    using Normal_3 = GEO::vec3;
    using Surface_patch_index = unsigned;
    unsigned nb_faces() const { return mesh.facets.nb(); }

    GEO::index_range face_range() const { return GEO::index_range{mesh.facets.begin(), mesh.facets.end()}; }

    unsigned patch_id(Face_descriptor f) const { return f; }
    unsigned nb_face_vertices(Face_descriptor face) const { return mesh.facets.nb_vertices(face); }
    auto face_vertices(Face_descriptor face) const { return mesh.facets.vertices(face); }
public:
    GEO::Mesh const &mesh;
};

class PolyLine_wrapper {
public:
    using Edge_descriptor = int;
    using Curve_index = unsigned;
    unsigned nb_edges() const { return mesh.edges.nb(); }
    GEO::index_range edge_range() const {
        return GEO::index_range{mesh.edges.begin(), mesh.edges.end()};
    }
    unsigned curve_id(Edge_descriptor e) const { return e; }
    int edge_vertex(Edge_descriptor edge, unsigned i) const { return mesh.edges.vertex(edge, i); }
public:
    GEO::Mesh const &mesh;
};


int main(int argc, char** argv) {
    GEO::initialize();
    const std::string filename = (argc > 1) ? argv[1] : "../data/fandisk_kenshi_hexmesh.mesh";
    const std::string surface_filename = (argc > 2) ? argv[2] : "";
    const std::string curves_filename = (argc > 3) ? argv[3] : "";
    const std::string point_targets_filename = (argc > 4) ? argv[4] : "";

    GEO::Mesh mesh;
    if(!GEO::mesh_load(filename, mesh)) {
        std::cerr << "Error loading mesh: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    GEO::Mesh reference;
    reference.copy(mesh);

    GEO::Mesh surface;
    if (!surface_filename.empty()) {
        if(!GEO::mesh_load(surface_filename, surface)) {
            std::cerr << "Error loading surface mesh: " << surface_filename << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        surface.copy(mesh);
    }
    surface.edges.clear();
    surface.cells.clear();
    surface.vertices.remove_isolated();
    surface.facets.triangulate();
    GEO::MeshFacetsAABB aabb;
    aabb.initialize(surface);

    std::vector<GEO::vec3> normal(surface.facets.nb());
    for (auto f : surface.facets) {
        normal[f] = GEO::Geom::mesh_facet_normal(surface, f);
        if (normal[f].length() > 1e-10)
            normal[f] /= normal[f].length();
        else
            normal[f] = GEO::vec3(0., 0., 1.);
    }

    auto query = [&](GEO::vec3 const &coord, unsigned surface_id, double radius) -> std::tuple<GEO::vec3, GEO::vec3, double> {
        GEO::vec3 res;
        double sqr_dist;
        GEO::index_t facet = aabb.nearest_facet(coord, res, sqr_dist);
        return {res, normal[facet], 1.};
    };


    GEO::Mesh curves;
    if (!curves_filename.empty()) {
        if(!GEO::mesh_load(curves_filename, curves)) {
            std::cerr << "Error loading curves mesh: " << curves_filename << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        curves.copy(mesh);
    }
    curves.facets.clear();
    curves.cells.clear();
    curves.vertices.remove_isolated();
    
    GEO::mesh_save(mesh, "input.mesh");
    GEO::mesh_save(surface, "surf.mesh");
    GEO::mesh_save(curves, "curves.mesh");

    std::vector<GEO::vec3> direction(curves.edges.nb());
    for (auto e : curves.edges) {
        curves.facets.create_triangle(curves.edges.vertex(e, 0), curves.edges.vertex(e, 1), curves.edges.vertex(e, 1));
        direction[e] = curves.vertices.point(curves.edges.vertex(e,1)) - curves.vertices.point(curves.edges.vertex(e,0));
        if (direction[e].length() > 1e-10) direction[e] /= direction[e].length();
        else direction[e] = GEO::vec3(0., 0., 1.);
    }
    GEO::MeshFacetsAABB aabb_curves;
    aabb_curves.initialize(curves);

    auto curve_query = [&](GEO::vec3 const &coord, unsigned edge_id, double radius) -> std::tuple<GEO::vec3, GEO::vec3, double> {
        GEO::vec3 res;
        double sqr_dist;
        GEO::index_t edge = aabb_curves.nearest_facet(coord, res, sqr_dist);
        return {res, direction[edge], 1.};
    };

    std::vector<std::tuple<GEO::index_t, GEO::vec3, double>> point_targets;

    if (!point_targets_filename.empty()) {
        std::ifstream infile(point_targets_filename);
        if (!infile) {
            std::cerr << "Error opening point targets file: " << point_targets_filename << std::endl;
            return EXIT_FAILURE;
        }
        GEO::index_t vertex_id;
        double x, y, z;
        while (infile >> vertex_id >> x >> y >> z) {
            point_targets.emplace_back(vertex_id, GEO::vec3(x, y, z), 1.);
        }
    }


    Mesh_wrapper mesh_wrapper(mesh, nullptr); // replacing nullptr by &reference will use it as a reference;
    mesh_wrapper.set_orientation(false, true, false, true); // beware of the orientation of your input elements!

    Boundary_wrapper surface_wrapper{mesh};
    PolyLine_wrapper curve_wrapper{mesh};

    Mesh_smoothing_3::Mesh_smoother smoother(mesh_wrapper, surface_wrapper, curve_wrapper);

    smoother.set_boundary_query(query);
    smoother.set_curves_query(curve_query);
    smoother.set_vertex_target_positions(point_targets);

    smoother.set_boundary_weight(Mesh_smoothing_3::Parameters::STRONG);
    smoother.set_verbose();
    smoother.set_max_number_of_iteration(100);
    smoother.run();

    GEO::mesh_save(mesh, "output.mesh");

    return EXIT_SUCCESS;
}
