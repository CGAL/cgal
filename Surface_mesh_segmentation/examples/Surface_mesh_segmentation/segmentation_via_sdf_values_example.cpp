#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input("data/cactus.off");
    if ( !input || !(input >> mesh) || mesh.empty() ) {
        std::cerr << "Not a valid off file." << std::endl;
        return EXIT_FAILURE;
    }

    // create a property-map for segment-ids
    typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

    // calculate SDF values and segment the mesh using default parameters.
    std::size_t number_of_segments = CGAL::segmentation_via_sdf_values(mesh, segment_property_map);

    std::cout << "Number of segments: " << number_of_segments << std::endl;

    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
        std::cout << segment_property_map[facet_it] << " ";
    }
    std::cout << std::endl;
}
