#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

int main()
{
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input("data/cactus.off");
    if ( !input || !(input >> mesh) || mesh.empty() || ( !CGAL::is_triangle_mesh(mesh)) ) {
        std::cerr << "Input is not a triangle mesh" << std::endl;
        return EXIT_FAILURE;
    }

    // create a property-map for segment-ids
    typedef std::map<face_descriptor, std::size_t> Face_int_map;
    Face_int_map internal_segment_map;
    boost::associative_property_map<Face_int_map> segment_property_map(internal_segment_map);

    // calculate SDF values and segment the mesh using default parameters.
    std::size_t number_of_segments = CGAL::segmentation_via_sdf_values(mesh, segment_property_map);

    std::cout << "Number of segments: " << number_of_segments << std::endl;

    // print segment-ids
    for(face_descriptor f : faces(mesh) ) {
        std::cout << segment_property_map[f] << " ";
    }
    std::cout << std::endl;
    return EXIT_SUCCESS;
}
