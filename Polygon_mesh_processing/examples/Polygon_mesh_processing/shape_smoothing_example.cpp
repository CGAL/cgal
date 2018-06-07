#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int main(int argc, char* argv[]){

    namespace PMP = CGAL::Polygon_mesh_processing;
    const char* filename = argc > 1 ? argv[1] : "data/pig.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    const double time = 0.1;
    const unsigned int nb_iterations = 2;

    PMP::smooth_along_curvature_flow(mesh, time,
        PMP::parameters::number_of_iterations(nb_iterations));

    std::ofstream output("mesh_shape_smoothed.off");
    output << mesh;
    output.close();

    return 0;
}
