#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int main(){

    namespace PMP = CGAL::Polygon_mesh_processing;
    const char* filename = "data/pig.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    /*
    unsigned int nb_iter = 10;
    CGAL::Polygon_mesh_processing::smooth_curvature_flow(mesh,
                                                         CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));

    std::ofstream output("data/dino_smoothed.off");
    output << mesh;
    output.close();

    */
    unsigned int nb_iter2 = 1;
    const double time = 10;
    CGAL::Polygon_mesh_processing::smooth_modified_curvature_flow(faces(mesh), mesh,
                                                                  CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter2));

    std::ofstream output2("data/pig_smoothed_to_sphere.off");
    output2 << mesh;
    output2.close();

    return 0;
}
