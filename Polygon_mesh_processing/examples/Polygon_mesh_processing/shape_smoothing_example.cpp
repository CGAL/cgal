#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>


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

    const double time = 1e-5;
    CGAL::Polygon_mesh_processing::smooth_modified_curvature_flow(mesh, time);

    std::ofstream output("data/pig_smoothed_to_sphere.off");
    output << mesh;
    output.close();

    return 0;
}
