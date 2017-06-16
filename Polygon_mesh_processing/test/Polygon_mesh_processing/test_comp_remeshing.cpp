#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/angle_smoothing.h>
#include <CGAL/Polygon_mesh_processing/area_smoothing.h>


#define CGAL_TEST_COMP_REMESHING_DEBUG
//#define CGAL_TEST_COMP_REMESHING_OUTPUT

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;



int main(int argc, char* argv[]){



    const char* filename;
    std::ifstream input;
    Mesh mesh;
#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    std::ofstream output;
#endif

    ///
    filename = "data/polygon3D.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());
    //CGAL::Polygon_mesh_processing::area_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default(), faces(mesh));


    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/polygon3D_smoothed.off");
    output << mesh;
    output.close();
#endif


    ///
    filename = "data/blobby_3cc.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());
    //CGAL::Polygon_mesh_processing::area_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default(), faces(mesh));


    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/blobby_3cc_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/cube_quad.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());
    //CGAL::Polygon_mesh_processing::area_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default(), faces(mesh));


    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/cube_quad_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/elephant.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/elephant_smoothed.off");
    output << mesh;
    output.close();
#endif

    filename = "data/degenerate_polygon.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/degenerate_polygon_smoothed.off");
    output << mesh;
    output.close();
#endif

    filename = "data/sneaky_degenerate_polygon.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/sneaky_degenerate_polygon_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/joint_refined.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/joint_refined_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/mannequin-devil.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();
#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/mannequin-devil_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/mech-holes-shark.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/mech-holes-shark_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/non_manifold_vertex.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }
    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/non_manifold_vertex_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/overlapping_triangles.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/overlapping_triangles_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/tetra1.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/tetra1_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/tetra2.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/tetra2_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/tetra3.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/tetra3_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/tetra4.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/tetra4_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/two_tris_collinear.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/two_tris_collinear_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/degtri_nullface.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/degtri_nullface_smoothed.off");
    output << mesh;
    output.close();
#endif

    ///
    filename = "data/U.off";
    input.open(filename);

#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }

    CGAL::Polygon_mesh_processing::angle_remeshing(mesh, CGAL::Polygon_mesh_processing::parameters::all_default());

    input.close();

#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open("data/U_smoothed.off");
    output << mesh;
    output.close();
#endif


    return 0;
}
