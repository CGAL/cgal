#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>


#define CGAL_TEST_COMP_REMESHING_DEBUG
#define CGAL_TEST_COMP_REMESHING_OUTPUT


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;



int main(int argc, char* argv[]){



    std::string filename;
    std::ifstream input;
    Mesh mesh;
#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    std::ofstream output;
#endif



    std::vector<std::string> filenames = {
        "data/polygon",
        "data/polygon3D",
        "data/blobby_3cc",
        "data/cube_quad",
        "data/elephant",
        "data/degenerate_polygon",
        "data/sneaky_degenerate_polygon",
        "data/joint_refined",
        "data/mannequin-devil",
        "data/mech-holes-shark",
        "data/non_manifold_vertex",
        "data/overlapping_triangles",
        "data/tetra1",
        "data/tetra2",
        "data/tetra3",
        "data/tetra4",
        "data/two_tris_collinear",
        "data/U"
    };


    for(auto i=0; i!= filenames.size(); ++i)
    {
        filename = filenames[i]+".off";
        input.open(filename);


#ifdef CGAL_TEST_COMP_REMESHING_DEBUG
std::cout<<"case: "<< filename << std::endl;
#endif

        if (!input || !(input >> mesh) || mesh.is_empty()) {
            std::cerr << "Not a valid .off file." << std::endl;
            return 1;
        }
        input.close();


        CGAL::Polygon_mesh_processing::angle_remeshing(mesh, faces(mesh), edges(mesh), CGAL::Polygon_mesh_processing::parameters::all_default());
        CGAL::Polygon_mesh_processing::area_remeshing(mesh, faces(mesh), edges(mesh), CGAL::Polygon_mesh_processing::parameters::all_default());



#ifdef CGAL_TEST_COMP_REMESHING_OUTPUT
    output.open(filenames[i]+"_smoothed"+".off");
    output << mesh;
    output.close();
#endif

    }



    return 0;
}
