#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char **argv)
{
    if (argc !=2){
        std::cerr << "Usage: " << argv[0] << " input.OFF" << std::endl;
        return 1;
    }
  
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input(argv[1]);
    
    if ( !input || !(input >> mesh) || mesh.empty() ){
        std::cerr << argv[1] << " is not a valid off file." << std::endl;
        return 1;
    }

    // create a property-map for segment-ids (it is an adaptor for this case)
    typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

    // calculate SDF values and segment the mesh using default parameters.
    int number_of_segments = CGAL::compute_sdf_values_and_segment(mesh, segment_property_map);

    std::cout << "Number of segments: " << number_of_segments << std::endl;

    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
        std::cout << segment_property_map[facet_it] << std::endl;                                 
    }
}