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
    const char* filename = "data/dino3.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }


    // SHAPE SMOOTHING
    double time = 1e-5;
    unsigned int nb_iter = 1;
    CGAL::Polygon_mesh_processing::smooth_shape(mesh, time, nb_iter);


    // MESH SMOOTHING
    CGAL::Polygon_mesh_processing::angle_smoothing(mesh);
    CGAL::Polygon_mesh_processing::area_smoothing(mesh);



    //CGAL::Polygon_mesh_processing::angles_evaluation<Mesh, K>(mesh, "data/angles_after.dat");
    //CGAL::Polygon_mesh_processing::areas_evaluation<Mesh, K>(mesh, "data/areas_after.dat");
    //CGAL::Polygon_mesh_processing::aspect_ratio_evaluation(mesh, "data/aspect_ratios_after.dat");

    std::ofstream output("data/dino3_smoothed.off");
    output << mesh;
    output.close();

    std::cout<<"OK."<<std::endl;


    return 0;
}
