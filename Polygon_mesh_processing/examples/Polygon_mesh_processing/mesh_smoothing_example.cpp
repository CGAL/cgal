#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>
#define CGAL_PMP_SMOOTHING_VERBOSE

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int main(int argc, char* argv[]){

    namespace PMP = CGAL::Polygon_mesh_processing;
    const char* filename = argc > 1 ? argv[1] : "data/sphere.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
      std::cerr << "Not a valid .off file." << std::endl;
      return 1;
    }

    unsigned int nb_iterations = 10;

    for(unsigned int t = 0 ; t < nb_iterations; ++t)
    {
      CGAL::Polygon_mesh_processing::smooth_angles(mesh);
      CGAL::Polygon_mesh_processing::smooth_areas(mesh);
      CGAL::Polygon_mesh_processing::smooth_angles(mesh);
    }

    std::ofstream output("data/sphere_smoothed.off");
    output << mesh;
    output.close();


    return 0;
}
