#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int main(int argc, char* argv[]){

    const char* filename = "data/mannequin-devil.off";
    std::ifstream input(filename);
    std::ofstream output;

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::curvature_flow(mesh);

    output.open("data/mannequin-devil_smoothed.off");
    output << mesh;
    output.close();


    return 0;
}
