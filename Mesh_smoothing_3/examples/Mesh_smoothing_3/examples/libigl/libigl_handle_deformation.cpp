#include <cstdlib>
#include <random>

#include <Mesh_smoothing_3/Mesh_smoothing_3.h>

#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include <iostream>

using Mesh_smoothing_3::utils::Contiguous_unsigned_range;

class Tetrahedral_mesh_wrapper {
public:
    using Cell_descriptor = unsigned;
    using Vertex_descriptor = int;
    using Point_3 = Eigen::Vector3d;

    unsigned nb_cells() const { return T.rows(); }
    unsigned nb_vertices() const { return V.rows(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return V.row(vertex); }
    void set_new_vertex_coordinates(int vertex, Point_3 coord) { V.row(vertex) = coord; }

    Contiguous_unsigned_range cell_range() const {
        return Contiguous_unsigned_range{0, nb_cells()};
    }
    std::array<int, 4> cell_vertices(Cell_descriptor cell) const {
        return {T(cell, 0), T(cell,1), T(cell,2), T(cell,3)};
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor cell) const {
        return {ref_V.row(T(cell, 0)), ref_V.row(T(cell,1)), ref_V.row(T(cell,2)), ref_V.row(T(cell,3))};
    }
public:
    Eigen::MatrixXd &V;
    Eigen::MatrixXd const &ref_V;
    Eigen::MatrixXi const &T;
};

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.mesh input_handles.txt" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string filename = argv[1];
    const std::string handlename = argv[2];

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXi F;

    igl::readMESH(filename, V, T, F);

    // Print the vertices and faces matrices
    std::cout << "Vertices: " << std::endl << V.rows() << std::endl;
    std::cout << "Tets:    " << std::endl << T.rows() << std::endl;
    std::cout << "Faces:    " << std::endl << F.rows() << std::endl;

    std::vector<std::pair<unsigned, Eigen::Vector3d>> handles;
    {
        std::ifstream handle_file(handlename);
        if (!handle_file.is_open()) {
            std::cerr << "Error opening handle file " << handlename << std::endl;
            return EXIT_FAILURE;
        }
        unsigned vid;
        double x, y, z;
        while (handle_file >> vid >> x >> y >> z) {
            handles.push_back({vid, Eigen::Vector3d{x, y, z}});
        }
    }

    igl::writeMESH("handle_input.mesh", V, T, F);

    Eigen::MatrixXd ref_V = V;


    Tetrahedral_mesh_wrapper mesh_wrapper {V, ref_V, T};

    Mesh_smoothing_3::Mesh_smoother smoother(mesh_wrapper);
    smoother.set_verbose();
    smoother.set_locked_boundary(false);

    smoother.set_vertex_target_positions(handles);

    bool res = smoother.run();

    std::cerr << "Smoothing result = " << (res ? "SUCCESS" : "FAILURE") << std::endl;
    igl::writeMESH("handle_output.mesh", V, T, F);
}