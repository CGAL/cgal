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
    const char* filename = "data/dino-coarse.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angles_evaluation<Mesh, K>(mesh, "data/angles_before.dat");
    CGAL::Polygon_mesh_processing::areas_evaluation<Mesh, K>(mesh, "data/areas_before.dat");
    //CGAL::Polygon_mesh_processing::aspect_ratio_evaluation(mesh, "data/aspect_ratios_before.dat");

    //CGAL::Polygon_mesh_processing::compatible_smoothing(mesh, PMP::parameters::protect_constraints(true));

    double time = 10 * 1e-6;
    //CGAL::Polygon_mesh_processing::smooth_shape(mesh, time, 1);

    CGAL::Polygon_mesh_processing::area_smoothing(mesh);

    CGAL::Polygon_mesh_processing::angles_evaluation<Mesh, K>(mesh, "data/angles_after.dat");
    CGAL::Polygon_mesh_processing::areas_evaluation<Mesh, K>(mesh, "data/areas_after.dat");
    //CGAL::Polygon_mesh_processing::aspect_ratio_evaluation(mesh, "data/aspect_ratios_after.dat");

    std::ofstream output("data/dino_smoothed.off");
    output << mesh;
    output.close();


    return 0;
}
