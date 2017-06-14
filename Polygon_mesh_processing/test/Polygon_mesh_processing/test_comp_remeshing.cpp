#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/angle_smoothing.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;



int main(int argc, char* argv[]){


    std::cout<<"hello\n";

    const char* filename;
    std::ifstream input;
    std::ofstream output;
    Mesh mesh;


    ///
    filename = "data/blobby_3cc.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_blobby_3cc.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/cube_quad.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_cube_quad.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/elephant.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_elephant.off");
    //output << mesh;
    //output.close();


    filename = "data/degenerate_polygon.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_degenerate_polygon.off");
    //output << mesh;
    //output.close();


    filename = "data/sneaky_degenerate_polygon.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_sneaky_degenerate_polygon.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/joint_refined.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_joint_refined.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/mannequin-devil.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_mannequin-devil.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/mech-holes-shark.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_mech-holes-shark.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/non_manifold_vertex.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }
    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_non_manifold_vertex.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/overlapping_triangles.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_overlapping_triangles.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/tetra1.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_tetra1.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/tetra2.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_tetra2.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/tetra3.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_tetra3.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/tetra4.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_tetra4.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/two_tris_collinear.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_two_tris_collinear.off");
    //output << mesh;
    //output.close();


    ///
    filename = "data/U.off";
    input.open(filename);


    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

    //output.open("data/smoothed_U.off");
    //output << mesh;
    //output.close();




    return 0;
}
