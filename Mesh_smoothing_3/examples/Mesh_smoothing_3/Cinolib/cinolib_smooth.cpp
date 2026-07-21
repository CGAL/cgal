#include <Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/octree.h>

using Mesh_smoothing_3::utils::Contiguous_unsigned_range;

class Mesh_wrapper : public Mesh_smoothing_3::helper_structures::Mixed_element_mesh<int, int, cinolib::vec3d, Contiguous_unsigned_range> {
public:
    std::size_t nb_vertices() const override { return mesh.num_verts(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const override {
        return mesh.vert(vertex);
    }
    void set_new_vertex_coordinates(int vertex, Point_3 coord) override { mesh.vert(vertex) = coord; }

    Contiguous_unsigned_range input_cell_range() const override { return Contiguous_unsigned_range{0, mesh.num_polys()}; }

    Vertex_descriptor get_cell_vertex(int cell, unsigned local_Vertex_descriptor) const override {
        return mesh.poly_vert_id(cell,local_Vertex_descriptor);
    };

    cinolib::vec3d get_ref_vertex_coordinates(int vertex) const override {
        assert(ref_mesh != nullptr);
        return ref_mesh->vert(vertex);
    }

    Shape const * get_element_shape(int cell) const override {
        switch (mesh.adj_p2v(cell).size()) {
            case 4:
                return &tet_ref;
                break;
            case 8:
                return &hex_ref;
                break;
            case 5:
                return &py_ref;
                break;
            case 6:
                return &we_ref;
                break;
            default:
                return nullptr;
                break;
        }
    }

public:
    Mesh_wrapper(cinolib::Polyhedralmesh<> &mesh_, cinolib::Polyhedralmesh<> const * ref_mesh_)
    : mesh(mesh_)
    , ref_mesh(ref_mesh_)
    {
        this->has_reference_mesh = (ref_mesh != nullptr);
        this->assemble(); // CRITICAL
    }

    void set_orientation(bool inv_tet, bool inv_hex, bool inv_pyr, bool inv_wed) {
        tet_ref.inverse = inv_tet;
        hex_ref.inverse = inv_hex;
        py_ref.inverse  = inv_pyr;
        we_ref.inverse  = inv_wed;
    }

    cinolib::Polyhedralmesh<> &mesh;
    cinolib::Polyhedralmesh<> const * ref_mesh;

    Mesh_smoothing_3::Shapes::VTK_TETRAHEDRON<cinolib::vec3d> tet_ref;
    Mesh_smoothing_3::Shapes::VTK_HEXAHEDRON<cinolib::vec3d> hex_ref;
    Mesh_smoothing_3::Shapes::VTK_PYRAMID<cinolib::vec3d> py_ref;
    Mesh_smoothing_3::Shapes::VTK_WEDGE<cinolib::vec3d> we_ref;
};

void compute_octree_from_input_surface(cinolib::Trimesh<> const &target_mesh, cinolib::Octree &octree, std::vector<cinolib::vec3d> &face_normal) {
    octree.build_from_mesh_polys(target_mesh);
    for(uint fid=0; fid<target_mesh.num_polys(); ++fid) {
        auto face_vertices = target_mesh.poly_verts_id(fid);
        cinolib::vec3d p0 = target_mesh.vert(face_vertices[0]);
        cinolib::vec3d p1 = target_mesh.vert(face_vertices[1]);
        cinolib::vec3d p2 = target_mesh.vert(face_vertices[2]);
        cinolib::vec3d n = (p1 - p0).cross(p2 - p0);
        n.normalize();
        face_normal.push_back(n);
    }
}

void compute_octree_from_mesh_surface(cinolib::Polyhedralmesh<> const &mesh, cinolib::Octree &octree, std::vector<cinolib::vec3d> &face_normal) {
    // extracting the surface for projection
    for(uint fid=0; fid<mesh.num_faces(); ++fid) {
        if (mesh.face_is_on_srf(fid)) {
            auto face_vertices = mesh.face_verts_id(fid);
            if (face_vertices.size() > 3) {
                cinolib::vec3d center = mesh.face_centroid(fid);
                for (uint i = 0; i < face_vertices.size(); ++i) {
                    cinolib::vec3d p0 = mesh.vert(face_vertices[(i+0)%face_vertices.size()]);
                    cinolib::vec3d p1 = mesh.vert(face_vertices[(i+1)%face_vertices.size()]);
                    cinolib::vec3d n = (p1-p0).cross(center-p0);
                    n.normalize();
                    octree.push_triangle(face_normal.size(), p0, p1, center);
                    face_normal.push_back(n);
                }
            }
            else if (face_vertices.size() == 3) {
                cinolib::vec3d p0 = mesh.vert(face_vertices[0]);
                cinolib::vec3d p1 = mesh.vert(face_vertices[1]);
                cinolib::vec3d p2 = mesh.vert(face_vertices[2]);
                cinolib::vec3d n = (p1 - p0).cross(p2 - p0);
                n.normalize();
                octree.push_triangle(face_normal.size(), p0, p1, p2);
                face_normal.push_back(n);
            }
            else {
                assert(false && "Unsupported face size");
            }
        }
    }
    octree.build();
}

int main(int argc, char** argv) {
    const std::string filename = (argc > 1) ? argv[1] : "../data/fandisk_kenshi_hexmesh.mesh";

    cinolib::Polyhedralmesh<> mesh(filename.c_str());
    Mesh_smoothing_3::helper_structures::Polygonal_boundary<unsigned, unsigned, cinolib::vec3d> boundary;
    for(uint fid=0; fid<mesh.num_faces(); ++fid) {
        if (mesh.face_is_on_srf(fid)) {
            auto face_vertices = mesh.face_verts_id(fid);
            boundary.add_polygon(face_vertices);
        }
    }

    cinolib::Octree octree;
    std::vector<cinolib::vec3d> face_normal;

    mesh.save("input.mesh");


    if (argc > 2) {
        cinolib::Trimesh<> target_mesh(argv[2]);
        compute_octree_from_input_surface(target_mesh, octree, face_normal);
    }
    else
    {
        compute_octree_from_mesh_surface(mesh, octree, face_normal);
    }

    auto query = [&](cinolib::vec3d pt, unsigned, double) -> std::tuple<cinolib::vec3d, cinolib::vec3d, double> {
        unsigned id;
        cinolib::vec3d proj;
        double dist;
        octree.closest_point(pt, id, proj, dist);
        return {proj, face_normal[id], 1.};
    };


    Mesh_wrapper mesh_wrapper(mesh, nullptr); // replacing nullptr by &reference will use it as a reference;
    mesh_wrapper.set_orientation(false, false, false, false); // beware of the orientation of your input elements!

    Mesh_smoothing_3::Mesh_smoother smoother(mesh_wrapper, boundary);

    smoother.set_boundary_query(query);

    smoother.set_verbose();
    smoother.set_max_number_of_iteration(100);
    smoother.run();

    mesh.save("output.mesh");


    return EXIT_SUCCESS;
}
