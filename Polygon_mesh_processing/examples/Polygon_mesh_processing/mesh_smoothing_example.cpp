#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>
#define CGAL_PMP_SMOOTHING_VERBOSE

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int main(){

    namespace PMP = CGAL::Polygon_mesh_processing;
    const char* filename = "data/sphere-api.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    unsigned int nb_iter = 10;


    for(int t = 0 ; t< 10; ++t)
    {
        CGAL::Polygon_mesh_processing::smooth_angles(mesh);
        CGAL::Polygon_mesh_processing::smooth_areas(mesh);
        CGAL::Polygon_mesh_processing::smooth_angles(mesh);

    }


    /*
    CGAL::Polygon_mesh_processing::compatible_smoothing(mesh,
                                                        CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));
                                                        */

    std::ofstream output("data/sphere-api_smoothed.off");
    output << mesh;
    output.close();


    return 0;
}
